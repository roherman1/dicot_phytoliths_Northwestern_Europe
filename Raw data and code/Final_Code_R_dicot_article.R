# 1. Load packages and data ----
library(readxl)    # For reading Excel files
library(dplyr)     # Data manipulation
library(tibble)    # Modern data frames
library(xlsx)      # Reading and writing Excel files
library(vegan)     # Community ecology statistics
library(lme4)      # Linear mixed-effects models
library(report)    # Automated reports of statistical models
library(ggplot2)   # Data visualisation
library(MASS)      # Statistical functions
library(ggforce)   # Enhanced ggplot2 visualisations
library(ggrepel)   # Repel overlapping text labels
library(ggpubr)    # Publication ready plots

# Load datasets
sample_table <- readRDS("sample_table.rdata")
phytolith_counts_table <- readRDS("phytolith_counts_table.rdata")

# 2. Calculate relative frequencies (RF) (%) ----
# * 2.1 Data pre-processing ----
# select rows
select_rows <- subset(phytolith_counts_table,
                      Info == "diagnostic"
                      |Info == "non-diagnostic"
                      |Info == "Name")

# Set the first row as column names
colnames(select_rows) <- as.character(unlist(select_rows[1, , drop = FALSE]))

# Remove the first row from the data frame
select_rows <- select_rows[-1, , drop = FALSE]

# Split info columns and raw counts
info_columns <- c("Name", "Morphotypes", "Compound_variable")
columns <- select_rows[info_columns]
counts <- select_rows[, !(names(select_rows) %in% info_columns)]
# Raw counts: Change make data.frame numeric
counts <- as.data.frame(lapply(counts, as.numeric))

# Raw counts: Change NA with zeros 
counts[is.na(counts)] <- 0

# * 2.2 Calculate RF ----
relative_frequencies <- apply(counts,2, function(x) x/sum(x)*100)
relative_frequencies_morphotypes <- cbind(columns, relative_frequencies)

# * 2.3 Export RF table (%) ----
write.xlsx(relative_frequencies_morphotypes, file = "Relative_frequencies_in_percentages.xlsx")

# 3. GLMER ----
# * 3.1 preparing data  ----

# factorise production level
sample_table$NP.T.CA <- factor(sample_table$NP.T.CA, levels = c("Non-producer", "Trace", "Common/abundant"), labels = c("Non-producer", "Trace", "Common/abundant"))
sample_table$NP.PRODUCER <- factor(sample_table$NP.PRODUCER, levels = c("NP", "PRODUCER"), labels = c("Non-producer", "Producer"))

# * 3.2 modelling ----
model_glmer <-glmer (NP.PRODUCER ~ Superorder + plant_part + (1|Specimen),
                     data = sample_table,
                     family = binomial,
                     control = glmerControl(optimizer="bobyqa",
                                            optCtrl=list(maxfun=2e5))) #https://svmiller.com/blog/2018/06/mixed-effects-models-optimizer-checks/

# * 3.3 report ----
report(model_glmer)
print()

library(broom.mixed)
tidy_output <- tidy(model_glmer)
print(tidy_output)
write.xlsx(tidy_output, "glmer_output.xlsx")

# 4. Production index (PI) ----
# * 4.1 PI and plant parts ----

# Proportion chart (relative frequencies)
table(sample_table$NP.T.CA,sample_table$plant_part)

ggplot(sample_table, aes(x = NP.T.CA, y = ..count../sum(..count..), fill = plant_part)) +
  geom_bar() +
  labs(y = "Frequency", x= "Production Index", fill= "Plant part") +
  scale_fill_discrete(labels=c('Rep. part (fruits or seeds)','Rep. part (inflorescences)','Leaves', 'Stem'))

# * 4.2 PI and taxonomy at family level  ----
# reverse order of Family to plot alphabetically
sample_table$Family <- factor(sample_table$Family)
sample_table$Family <- factor(sample_table$Family, levels = rev(levels(sample_table$Family)))

#plot
ggplot(data = sample_table, aes(x = Family, fill = NP.T.CA)) +
  geom_bar(position = "fill") + 
  ylab("Proportion") + 
  theme_classic() +
  stat_count(geom = "text", 
             aes(label = stat(count)),
             position = position_fill(vjust = 0.5), 
             colour = "white", size = 3) +
  coord_flip() + 
  labs(fill = "Production levels")

# * 4.3 PI and life form ----

# Proportion chart (relative frequencies)
table(sample_table$NP.T.CA,sample_table$Plant_category, sample_table$plant_part)

ggplot(sample_table, aes(x = NP.T.CA, y = ..count../sum(..count..), fill = Plant_category )) +
  geom_bar() +
  labs(y = "Frequency", x= "Production Index", fill= "Life form") 

#5. Silicification ----

#* 5.1 Silicification and PI  ----

# select producing phytolith samples
silicification <- subset(sample_table, silicification == "well-silicified"|silicification == "weakly silicified")

plot1 <- ggplot(silicification, aes(x = silicification, fill = NP.T.CA)) +
  geom_bar(aes(y = ..count../sum(..count..)), width = 0.75) +
  theme_classic() + 
  labs(y = "Relative Frequency", x = "Silicification rate", fill = "Production Index") +
  theme(legend.position = "right")

#* 5.2 Silicification and taxonomy ----
plot2 <- ggplot(silicification, aes(x = silicification, fill = Superorder)) +
  geom_bar(aes(y = ..count../sum(..count..)), width = 0.75) +
  theme_classic() + 
  labs(y = "Relative Frequency", x = "Silicification rate", fill = "Superorder") +
  theme(legend.position = "right")

# reverse order of Family to plot alphabetically
sample_table$Family <- factor(sample_table$Family)
sample_table$Family <- factor(sample_table$Family, levels = rev(levels(sample_table$Family)))

# Set the levels in the desired order
sample_table$silicification <- ifelse(is.na(sample_table$silicification), "NA", sample_table$silicification)
sample_table$silicification <- factor(sample_table$silicification, 
                                      levels = c("NA", "weakly silicified", "well-silicified"))


# Plotting
ggplot(data = sample_table, aes(x = Family, fill = silicification)) +
  geom_bar(position = "fill") + 
  ylab("Proportion") + 
  theme_classic() +
  stat_count(geom = "text", 
             aes(label = stat(count)),
             position = position_fill(vjust = 0.5), 
             colour = "white", size = 3) +
  coord_flip() + 
  labs(fill = "Silicification rate")

#* 5.3 Silicification and plant part ----
plot3 <- ggplot(silicification, aes(x = silicification, fill = plant_part)) +
  geom_bar(aes(y = ..count../sum(..count..)), width = 0.75) +
  theme_classic() + 
  labs(y = "Relative Frequency", x = "Silicification rate", fill = "Plant part") +
  theme(legend.position = "right")

#* 5.3 Silicification and life form  ----
plot4 <- ggplot(silicification, aes(x = silicification, fill = Plant_category)) +
  geom_bar(aes(y = ..count../sum(..count..)), width = 0.75) +
  theme_classic() + 
  labs(y = "Relative Frequency", x = "Silicification rate", fill = "Plant part") +
  theme(legend.position = "right")

#* 5.4 Combined plots of PI, plant part, and taxonomy ----
arranged_plots <- ggarrange(plot1, plot2, plot3, plot4,
                            ncol = 2, nrow = 2,
                            common.legend = FALSE, legend = "right")

print(arranged_plots)

#* 5.5 Chi squared test ----
table(silicification$NP.T.CA, silicification$silicification)
table(silicification$Superorder, silicification$silicification)
table(silicification$plant_part, silicification$silicification)
table(silicification$Plant_category, silicification$silicification)

chi_squared_test1 <- chisq.test(silicification$NP.T.CA, silicification$silicification, correct = FALSE)
chi_squared_test2 <- chisq.test(silicification$Superorder, silicification$silicification, correct = FALSE)
chi_squared_test3 <- chisq.test(silicification$plant_part, silicification$silicification, correct = FALSE)
chi_squared_test4 <- chisq.test(silicification$Plant_category, silicification$silicification, correct = FALSE)

# Create an overview table
overview_table <- rbind(
  c("Production Index", chi_squared_test1$statistic, chi_squared_test1$p.value),
  c("Superorder", chi_squared_test2$statistic, chi_squared_test2$p.value),
  c("Plant Part", chi_squared_test3$statistic, chi_squared_test3$p.value),
  c("Life form", chi_squared_test4$statistic, chi_squared_test4$p.value))

# Print the overview table
colnames(overview_table) <- c("Variable", "Chi-squared", "p-value")
print(overview_table)
write.xlsx(overview_table, "chi_squared_results.xlsx")

# 6. Intra- and inter group variability in phytolith assemblages ----
# * 6.1 Data pre-processing ----
# select rows
select_rows_1 <- subset(phytolith_counts_table,
                      Info == "diagnostic"
                      |Info == 'Species_plant_part'
                      |Info =='200_diagnostic_phytoliths')

select_rows_1 <- t(select_rows_1)

# select samples with 200 diagnostic phytoliths
diagnostic_200_pre <- select_rows_1[select_rows_1[, 1] == "200 diagnostic phytoliths = YES", ]
diagnostic_200_pre <- diagnostic_200_pre[-(1:2), ]
labels <- diagnostic_200_pre[, 2]
diagnostic_200 <- diagnostic_200_pre[, -(1:2)]

# Convert the matrix to a data frame and make it numeric
diagnostic_200 <- as.data.frame(diagnostic_200)
diagnostic_200 <- as.data.frame(lapply(diagnostic_200, as.numeric))

# Change NA with zeros 
diagnostic_200[is.na(diagnostic_200)] <- 0

# * 6.2. Perform PERMANOVA  ----
# make a distance matrix
distance <-vegdist(diagnostic_200, method="horn") 
distance_matrix <- as.matrix(distance)
veg
# perform the permanova on the distance matrix

set.seed(1997)
adonis2 <- adonis2(distance_matrix ~ labels)
print(adonis2)
adonis_result_summary <- print(adonis2)
write.table(adonis_result_summary, file = "adonis_results.txt", quote = FALSE, sep = "\t")

write.table(labels, file = "labels.txt", quote = FALSE, sep = "\t")


# * 4.3 Dissimilarity intra-group data ----
# Calculate the average/sd dissimilarity
mean_value <- mean(distance_matrix)
print(mean_value)
sd_value <- sd(distance_matrix)
print(sd_value)

# prepare data frame with labels
labels_frame <-data.frame(labels)
diagnostic_200_labels<- cbind(labels, diagnostic_200)

# Identify instances that occur two times based on identical values in the first row
duplicate_instances <- diagnostic_200_labels[duplicated(diagnostic_200_labels$labels) | duplicated(diagnostic_200_labels$labels, fromLast = TRUE), ]

# Create an empty list to store Horn distances
distances_list <- list()

# Loop through each unique label that occurs two times
for (label in unique(duplicate_instances$labels)) {
  # Subset the data frame for the current label
  subset_data <- duplicate_instances[duplicate_instances$labels == label, -1]  # Exclude the 'labels' column
  
  # Calculate Horn's distance
  horn_distance <- vegdist(subset_data, method = "horn")
  
  # Store the result in the list
  distances_list[[label]] <- horn_distance
}

# Function to extract lower triangular part of a distance matrix
lower_triangular <- function(mat) {
  mat[lower.tri(mat, diag = TRUE)]
}

# Apply the function to each element in the list
lower_triangular_list <- lapply(distances_list, lower_triangular)
distances_df <- data.frame(do.call(cbind, lower_triangular_list))

# Clean the data.frame
pairwise_distances_duplicate_samples <- data.frame(t(distances_df[1,]))
rownames(pairwise_distances_duplicate_samples) <- gsub("\\.", " ", rownames(pairwise_distances_duplicate_samples))

pairwise_distances_duplicate_samples$RowNames <- rownames(pairwise_distances_duplicate_samples)
rownames(pairwise_distances_duplicate_samples) <- NULL

# Order the dataframe based on X1 in descending order
sorted_df <- pairwise_distances_duplicate_samples[order(pairwise_distances_duplicate_samples$X1), ]
sorted_df$RowNames <- factor(sorted_df$RowNames, levels = sorted_df$RowNames[order(sorted_df$X1)])

ggplot(sorted_df, aes(x = RowNames, y = X1)) +
  geom_bar(stat = "identity", color = "black", fill = "black") +
  geom_hline(yintercept = mean_value, linetype = "dotted", color = "red", size = 0.8, aes(linetype = "Mean")) +
  geom_hline(yintercept = mean_value + sd_value, linetype = "dotted", color = "blue", size = 0.8, aes(linetype = "SD")) +
  geom_hline(yintercept = mean_value - sd_value, linetype = "dotted", color = "blue", size = 0.8) +
  labs(title = "Pairwise dissimilarity values", x = "Duplicate samples with 200 diagnostic phytoliths", y = "Dissimilarity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  geom_text(aes(label = round(X1, 3)), vjust = -0.5, size = 2) +
  theme(legend.position = "bottom") +
  scale_linetype_manual(name = "Lines", values = c("Mean" = "dotted", "SD" = "dotted"), labels = c("Mean", "SD"))

# 7. Linear discriminant analysis ----
# * 7.1 Plant part differentiation ----
select_rows_2 <- subset(phytolith_counts_table,
                        Info == "diagnostic"
                        |Info == 'Plant part (seeds and fruits) together'
                        |Info =='200_diagnostic_phytoliths')

select_rows_2 <- t(select_rows_2)

# make labels_2: plant_part
diagnostic_200_pre_2 <- select_rows_2[select_rows_2[, 1] == "200 diagnostic phytoliths = YES", ]
diagnostic_200_pre_2 <- diagnostic_200_pre_2[-(1:2), ]
labels_2 <- diagnostic_200_pre_2[, 2]
labels_2 <- factor(labels_2)

# select samples with 200 diagnostic phytoliths
diagnostic_200_LDA <- diagnostic_200

select_rows_3 <- subset(phytolith_counts_table,
                        Info == "diagnostic")
result <- t(select_rows_3[, 2])

# remove columns with 0 or 1 morphotypes present 
colnames(diagnostic_200_LDA) <- result
LDA_data <- diagnostic_200[, colSums(diagnostic_200 != 0) > 1]
LDA_data_excluded <- diagnostic_200[, !(colnames(diagnostic_200) %in% colnames(LDA_data))]

# Save the original column names
original_colnames <- colnames(diagnostic_200_LDA)
# Rename the columns of LDA_data using the original column names
colnames(LDA_data) <- original_colnames[colSums(diagnostic_200 != 0) > 1]

# Tra-5 B, Tra-5 C, TRI_STE, and Cl-4 were removed (34 diagnostic morphotypes remaining) 

#scale the data
scaled_data_LDA <- as.data.frame(scale(LDA_data))
lda_model_1 <- lda(labels_2 ~ ., data = scaled_data_LDA)
lda_data_1 <- predict(lda_model_1, newdata = scaled_data_LDA)

# Create a data frame for plotting
plot_data_1 <- data.frame(
  LDA1 = lda_data_1$x[, 1],
  LDA2 = lda_data_1$x[, 2],
  Label = labels_2  # Assuming 'labels' is the column in your_data
)

# big plot
plot_1_plant_part <- ggplot(plot_data_1, aes(x = LDA1, y = LDA2, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label), type = "norm", level = 0.95, geom = "polygon", alpha = 0.2) +
  labs(title = "Grouping of plant part", x = "LDA1 (92%)", y = "LDA2 (6%)") +
  scale_color_discrete(name = "Group") +
  theme_minimal() + 
  annotate("rect", xmin = -4, xmax = 0, ymin = -6, ymax = 7, color = "black", fill = NA, linetype = 1)
  #geom_text(aes(label = labels), nudge_x = 0.1, nudge_y = 0.1, size = 3, check_overlap = TRUE) 

# sub plot
plot_1_plant_part_subplot <- ggplot(plot_data, aes(x = LDA1, y = LDA2, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label), type = "norm", level = 0.95, geom = "polygon", alpha = 0.2) +
  labs(title = "Subplot", x = "LDA1 (0.92%)", y = "LDA2 (0.06%)") +
  scale_color_discrete(name = "Group") +
  theme_minimal() +
  #geom_text(aes(label = labels), nudge_x = 0.1, nudge_y = 0.1, size = 3, check_overlap = TRUE) +
  coord_cartesian(xlim = c(-4, 0), ylim = c(-6, 7)) +
  annotate("rect", xmin = -4, xmax = 0, ymin = -6, ymax = 7, color = "black", fill = NA, linetype = 1)


# * 7.2 Life form differentation ----
select_rows_3 <- subset(phytolith_counts_table,
                        Info == "diagnostic"
                        |Info == 'Life form'
                        |Info =='200_diagnostic_phytoliths')

select_rows_3 <- t(select_rows_3)

# make labels_2: plant_part
diagnostic_200_pre_3 <- select_rows_3[select_rows_3[, 1] == "200 diagnostic phytoliths = YES", ]
diagnostic_200_pre_3 <- diagnostic_200_pre_3[-(1:2), ]
labels_3 <- diagnostic_200_pre_3[, 2]
labels_3 <- factor(labels_3)

# select samples with 200 diagnostic phytoliths
diagnostic_200_LDA <- diagnostic_200

select_rows_3 <- subset(phytolith_counts_table,
                        Info == "diagnostic")
result <- t(select_rows_3[, 2])

# remove columns with 0 or 1 morphotypes present 
colnames(diagnostic_200_LDA) <- result
LDA_data <- diagnostic_200[, colSums(diagnostic_200 != 0) > 1]
LDA_data_excluded <- diagnostic_200[, !(colnames(diagnostic_200) %in% colnames(LDA_data))]
# Tra-5 B, Tra-5 C, TRI_STE, and Cl-4 were removed (34 diagnostic morphotypes remaining) 

#scale the data
scaled_data_LDA <- as.data.frame(scale(LDA_data))
lda_model_2 <- lda(labels_3 ~ ., data = scaled_data_LDA)
lda_data_2 <- predict(lda_model_2, newdata = scaled_data_LDA)

# Create a data frame for plotting
plot_data_2 <- data.frame(
  LDA1 = lda_data_2$x[, 1],
  Label = labels_3,
  tags = labels
)

# Create a scatter plot with LD1 on the x-axis
plot_2_life_form <- ggplot(plot_data_2, aes(x = 1, y = LDA1, color = Label, label = labels)) +
  geom_point() +
  labs(title = "Grouping of life form", x = "", y = "LD1 (100%)") +
  scale_color_discrete(name = "Group") +
  theme_minimal() + 
  #geom_text_repel(box.padding = 0.5, segment.color = "grey", size = 2, max.overlaps = Inf) +
  theme(legend.position = "bottom") +
  annotate("rect", xmin = 0.98, xmax = 1.02, ymin = -0.3, ymax = 0.8, color = "black", fill = NA, linetype = 1) + 
  coord_cartesian(xlim = c(0.975, 1.025))


# Create the subplot with labels within coord_cartesian limits
plot_2_life_form_subplot <- subplot <- ggplot(plot_data_2, aes(x = 1, y = LDA1, color = Label, label = tags)) +
  geom_point() +
  labs(title = "Grouping of Life form", x = "", y = "LD1 (100%)") +
  scale_color_discrete(name = "Group") +
  theme_minimal() + 
  geom_text_repel(
    box.padding = 0.65,
    segment.color = "grey",
    size = 2,
    max.overlaps = Inf,
    data = subset(plot_data_2, LDA1 >= -0.3 & LDA1 <= 0.8)  # Filter labels within coord_cartesian limits
  ) +
  theme(legend.position = "bottom") +  # Move legend to the bottom
  coord_cartesian(xlim = c(0.975, 1.025), ylim = c(-0.6, 1.1)) +
  annotate("rect", xmin = 0.98, xmax = 1.02, ymin = -0.3, ymax = 0.8, color = "black", fill = NA, linetype = 1) 
  
# * 7.3 Taxonomic group differentation ----
select_rows_4 <- subset(phytolith_counts_table,
                        Info == "diagnostic"
                        |Info == 'Order'
                        |Info =='200_diagnostic_phytoliths')

select_rows_4 <- t(select_rows_4)
# make labels_2: plant_part
diagnostic_200_pre_4 <- select_rows_4[select_rows_4[, 1] == "200 diagnostic phytoliths = YES", ]
diagnostic_200_pre_4 <- diagnostic_200_pre_4[-(1:2), ]
labels_4 <- diagnostic_200_pre_4[, 2]
labels_4 <- factor(labels_4)
table(labels_4)
# select samples with 200 diagnostic phytoliths
diagnostic_200_LDA <- diagnostic_200

# remove columns with 0 or 1 morphotypes present 
colnames(diagnostic_200_LDA) <- result
LDA_data <- diagnostic_200[, colSums(diagnostic_200 != 0) > 1]
LDA_data_excluded <- diagnostic_200[, !(colnames(diagnostic_200) %in% colnames(LDA_data))]
# Tra-5 B, Tra-5 C, TRI_STE, and Cl-4 were removed (34 diagnostic morphotypes remaining) 

#scale the data
scaled_data_LDA <- as.data.frame(scale(LDA_data))
lda_model_3 <- lda(labels_4 ~ ., data = scaled_data_LDA)
lda_data_3 <- predict(lda_model_3, newdata = scaled_data_LDA)

# Create a data frame for plotting
plot_data_3 <- data.frame(
  LDA1 = lda_data_3$x[, 1],
  LDA2 = lda_data_3$x[, 2],
  Label = labels_4  # Assuming 'labels' is the column in your_data
)

# big plot
plot_3_taxonomy <- ggplot(plot_data_3, aes(x = LDA1, y = LDA2, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label), type = "norm", level = 0.95, geom = "polygon", alpha = 0.2) +
  labs(title = "Grouping of orders", x = "LDA1 (0.47%)", y = "LDA2 (0.30%)") +
  scale_color_discrete(name = "Group") +
  theme_minimal() +
  annotate("rect", xmin = -7, xmax = 7, ymin = -6, ymax = 5, color = "black", fill = NA, linetype = 1) 
  #geom_text(aes(label = labels), nudge_x = 5, nudge_y = 0.1, size = 1, check_overlap = TRUE) 

# sub plot
plot_3_taxonomy_subplot <- 
ggplot(plot_data_3, aes(x = LDA1, y = LDA2, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label), type = "norm", level = 0.95, geom = "polygon", alpha = 0.2) +
  labs(title = "Grouping of orders", x = "LDA1 (0.47%)", y = "LDA2 (0.30%)") +
  scale_color_discrete(name = "Group") +
  theme_minimal() +
  coord_cartesian(xlim = c(-7,7), ylim = c(-6,5)) +
  annotate("rect", xmin = -7, xmax = 7, ymin = -6, ymax = 5, color = "black", fill = NA, linetype = 1) 
  #geom_text(aes(label = labels), nudge_x = 5, nudge_y = 1, size = 3, check_overlap = TRUE) 

#* 7.4 comparing LDA results ----
# Extract eigenvalues from the three LDA models
eigenvalues_1 <- lda_model_1$svd^2
eigenvalues_2 <- lda_model_2$svd^2
eigenvalues_3 <- lda_model_3$svd^2

# Calculate the number of LDs
n_lds <- min(length(eigenvalues_1), length(eigenvalues_2), length(eigenvalues_3))

# Calculate the proportions of trace explained for all LDs
prop_explained_1 <- eigenvalues_1 / sum(eigenvalues_1)
prop_explained_2 <- eigenvalues_2 / sum(eigenvalues_2)
prop_explained_3 <- eigenvalues_3 / sum(eigenvalues_3)

# Combine all proportions into one list
prop_explained_list <- list(prop_explained_1, prop_explained_2, prop_explained_3)

# Find the maximum length in the list
max_length <- max(lengths(prop_explained_list))

# Pad shorter vectors with NA
prop_explained_padded <- lapply(prop_explained_list, function(x) c(x, rep(NA, max_length - length(x))))

# Create a data frame
prop_explained <- data.frame(do.call(cbind, prop_explained_padded))

# Print the data frame
print(prop_explained)

# Extract coefficients for LD1 from each model
coefficients_1 <- lda_model_1$scaling[, 1]
coefficients_2 <- lda_model_2$scaling[, 1]
coefficients_3 <- lda_model_3$scaling[, 1]

# Extract coefficients for LD2 from each model
coefficients_1_2 <- lda_model_1$scaling[, 2]
coefficients_3_2 <- lda_model_3$scaling[, 2]

# Create a data frame with coefficients
coefficients_table <- data.frame(
  Variable = colnames(scaled_data_LDA),  # Assuming scaled_data_LDA is the data used for LDA
  Plant_part_LD1 = coefficients_1,
  Plant_part_LD2 = coefficients_1_2,
  Life_form_LD1 = coefficients_2,
  Taxonomy_LD1 = coefficients_3,
  Taxonomy_LD2 = coefficients_3_2)



