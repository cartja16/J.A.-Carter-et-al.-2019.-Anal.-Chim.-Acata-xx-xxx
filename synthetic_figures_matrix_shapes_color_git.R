# Notes ====

# The following provides the code for PCA and AP Clustering to produce Figure 6

# Clear workspace
rm(list=ls())
dev.off()

# load packages
library(tidyverse)
library(readxl)
library(apcluster)
library(gridExtra)
library(ggrepel)
library(RColorBrewer)
library(viridis)
library(extrafont)

windowsFonts(Times=windowsFont("TT Times New Roman"))
par(family = "Times")

# Read in and format data ====

curve_points <- c("0 ppb", "100 ppb", "200 ppb", "500 ppb", "1000 ppb", "0 ppb_Sea", "100 ppb_Sea", "200 ppb_Sea", "500 ppb_Sea", "1000 ppb_Sea")

OES_Na_1_norm_wide <- read.csv(file = "G:/My Drive/WFU/BG_multivariate/R/BG_PCA/Na_1_norm.csv", header=TRUE, sep=",")
OES_Na_2_norm_wide <- read.csv(file = "G:/My Drive/WFU/BG_multivariate/R/BG_PCA/Na_2_norm.csv", header=TRUE, sep=",")
OES_Ca_1_norm_wide <- read.csv(file = "G:/My Drive/WFU/BG_multivariate/R/BG_PCA/Ca_1_norm.csv", header=TRUE, sep=",") %>%
  dplyr::filter(!Sample.Name %in% curve_points)
OES_Ca_2_norm_wide <- read.csv(file = "G:/My Drive/WFU/BG_multivariate/R/BG_PCA/Ca_2_norm.csv", header=TRUE, sep=",")
OES_Carb_1_norm_wide <- read.csv(file = "G:/My Drive/WFU/BG_multivariate/R/BG_PCA/Carb_1_norm.csv", header=TRUE, sep=",") %>%
  dplyr::filter(!Sample.Name %in% curve_points)
OES_Carb_2_norm_wide <- read.csv(file = "G:/My Drive/WFU/BG_multivariate/R/BG_PCA/Carb_2_norm.csv", header=TRUE, sep=",")
OES_Na_C_1_norm_wide <- read.csv(file = "G:/My Drive/WFU/BG_multivariate/R/BG_PCA/Na_C_1_norm.csv", header=TRUE, sep=",")

rid_carbon <- c("C.247.856", "C.193.027", "C.175.122", "Tb.308.958", "Ta.337.649", "N.174.465", "N.174.213")
analytes <- c("Cd.214.439", "Co.238.892", "Cr.267.716", "Pb.220.353")

matrices_all <- OES_Na_1_norm_wide %>%
  dplyr::bind_rows(OES_Na_2_norm_wide,
                   OES_Ca_1_norm_wide,
                   OES_Ca_2_norm_wide,
                   OES_Carb_1_norm_wide,
                   OES_Carb_2_norm_wide,
                   OES_Na_C_1_norm_wide) %>%
  dplyr::mutate(Sample.Name = case_when(Sample.Name == "0.1 % Na + TOC" ~ "0.1 % Na + 0.1 % TOC",
                                        Sample.Name == "0.2 % Na + TOC" ~ "0.2 % Na + 0.2 % TOC",
                                        Sample.Name == "0.1 % Ca 0.1 % TOC" ~ "0.1 % Ca + 0.1 % TOC",
                                        Sample.Name == "0.2 % Ca 0.2 % TOC" ~ "0.2 % Ca + 0.2 % TOC",
                                        !Sample.Name %in% c("0.1 % Na + TOC", "0.2 % Na + TOC", "0.1 % Ca 0.1 % TOC", "0.2 % Ca 0.2 % TOC") ~ Sample.Name))
matrices_all$Sample.Name <- gsub("%", "% m/v", matrices_all$Sample.Name)

Samples_levels <- c("0.1 % m/v TOC", "0.2 % m/v TOC", "1 % m/v TOC", "2 % m/v TOC",
                    "0.1 % m/v Na", "0.1 % m/v Na + 0.1 % m/v TOC", "0.2 % m/v Na", "0.2 % m/v Na + 0.2 % m/v TOC", "1 % m/v Na", "2 % m/v Na",
                    "0.1 % m/v Ca", "0.1 % m/v Ca + 0.1 % m/v TOC", "0.2 % m/v Ca", "0.2 % m/v Ca + 0.2 % m/v TOC")

Samples <- matrices_all %>%
  dplyr::mutate(Sample.Name = as.character(Sample.Name)) %>%
  dplyr::select(Sample.Name) %>%
  dplyr::filter(!Sample.Name %in% curve_points) %>%
  dplyr::pull()

matrices_all <- matrices_all %>%
  dplyr::mutate(Sample.Name = factor(Sample.Name, levels = c(curve_points, Samples_levels), ordered = TRUE))

BGVars <- matrices_all %>%
  dplyr::select(-Sample.Name, -analytes, - rid_carbon) %>%
  colnames()

# PCA stats ====

matrices_all <- matrices_all %>%
  dplyr::select(-analytes, -rid_carbon)

matrices_all_PCA <- matrices_all %>%
  dplyr::select(-Sample.Name) %>%
  scale(center = TRUE, scale = TRUE) %>%
  prcomp()

matrices_all_collect <- as.tibble(matrices_all_PCA$x) %>%
  dplyr::bind_cols(matrices_all)

PCA_avgs <- matrices_all_collect %>%
  dplyr::filter(Sample.Name %in% curve_points) %>%
  dplyr::summarize_at(vars(PC1, PC2), funs(mean(., na.rm=TRUE))) %>%
  as.numeric()

PC_Euclidean_tmp <- matrices_all_collect %>%
  dplyr::filter(Sample.Name %in% curve_points) %>%
  dplyr::mutate(Euclidean.Distance = "")

PC_Euclidean <- matrices_all_collect %>%
  dplyr::filter(!Sample.Name %in% curve_points) %>%
  dplyr::mutate(Euclidean.Distance = as.character(signif(sqrt((PC1-PCA_avgs[1])^2 + (PC2-PCA_avgs[2])^2), digits = 3))) %>%
  dplyr::bind_rows(PC_Euclidean_tmp)

# AP Clustering on PC scores ====

AP_cluster <- PC_Euclidean %>%
  dplyr::select(-BGVars, -Euclidean.Distance, -Sample.Name) 

AP_cluster_sim <- negDistMat(AP_cluster, r=2)
AP_cluster_model <- apcluster(AP_cluster_sim, details = TRUE)

AP_cluster_plot <- PC_Euclidean 

Exemplars_ind <- AP_cluster_model@exemplars

AP_cluster_plot$Exemplars <- "Non_Exemplar"

AP_cluster_plot$Exemplars[Exemplars_ind] <- "Exemplar"

AP_cluster_plot <- AP_cluster_plot %>%
  dplyr::mutate(Group = factor(AP_cluster_model@idx, levels = unique(AP_cluster_model@idx), ordered = TRUE)) %>%
  dplyr::mutate(Exemplars = factor(Exemplars, levels = c("Non_Exemplar", "Exemplar"))) %>%
  dplyr::mutate(Sample.Name = factor(Sample.Name, levels = c(curve_points, Samples_levels), ordered = TRUE)) %>%
  dplyr::group_by(Group) %>%
  dplyr::mutate(Exemp_x = PC1[Exemplars == "Exemplar"],
                Exemp_y = PC2[Exemplars == "Exemplar"]) %>%
  dplyr::ungroup()

# Read in and format recoveries ====

Na_matrix_recoveries <- read.csv(file = "G:/My Drive/WFU/BG_multivariate/R/BG_PCA/Na_matrix_recoveries.csv", header=TRUE, sep=",")
Ca_matrix_recoveries <- read.csv(file = "G:/My Drive/WFU/BG_multivariate/R/BG_PCA/Ca_matrix_recoveries.csv", header=TRUE, sep=",")
Na_C_matrix_recoveries <- read.csv(file = "G:/My Drive/WFU/BG_multivariate/R/BG_PCA/Na_C_matrix_recoveries.csv", header=TRUE, sep=",")
Carb_matrix_recoveries <- read.csv(file = "G:/My Drive/WFU/BG_multivariate/R/BG_PCA/Carb_matrix_recoveries.csv", header=TRUE, sep=",")
Sea_matrix_recoveries <- read.csv(file = "G:/My Drive/WFU/BG_multivariate/R/BG_PCA/Sea_recoveries.csv", header=TRUE, sep=",")

matrix_recoveries <- Na_matrix_recoveries %>%
  dplyr::bind_rows(Ca_matrix_recoveries,
                   Na_C_matrix_recoveries,
                   Carb_matrix_recoveries,
                   Sea_matrix_recoveries) %>%
  dplyr::mutate(Element = as.character(Element)) %>%
  dplyr::mutate(Percent_recoveries = (ugL_mean/500)*100) %>%
  dplyr::mutate(Matrix = case_when(Matrix == "0.1 % Na + TOC" ~ "0.1 % Na + 0.1 % TOC",
                                   Matrix == "0.2 % Na + TOC" ~ "0.2 % Na + 0.2 % TOC",
                                   Matrix == "0.1 % Ca 0.1 % TOC" ~ "0.1 % Ca + 0.1 % TOC",
                                   Matrix == "0.2 % Ca 0.2 % TOC" ~ "0.2 % Ca + 0.2 % TOC",
                                   !Matrix %in% c("0.1 % Na + TOC", "0.2 % Na + TOC", "0.1 % Ca 0.1 % TOC", "0.2 % Ca 0.2 % TOC") ~ Matrix)) %>%
  dplyr::mutate(RSD_percent = (ugL_sd/ugL_mean)*100)

matrix_recoveries$Matrix <- gsub("%", "% m/v", matrix_recoveries$Matrix)

recovery_range <- matrix_recoveries %>%
  dplyr::filter(Matrix %in% Samples) %>%
  dplyr::group_by(Matrix) %>%
  dplyr::summarize(min = min(Percent_recoveries),
                   max = max(Percent_recoveries),
                   mean = mean(Percent_recoveries)) %>%
  dplyr::mutate(mean_lab = paste(signif(mean, digits = 3), "%", sep = " "))

recovery_range$`Euclidean Distance` <- NA

for (i in seq_along(recovery_range$Matrix)) {
  recovery_range$`Euclidean Distance`[i] <- AP_cluster_plot$Euclidean.Distance[AP_cluster_plot$Sample.Name == recovery_range$Matrix[i]]
}

recovery_range$Matrix_label <- NA
recovery_range$mean_lab_let <- NA
recovery_range$dist.recov <- NA

recovery_range <- recovery_range %>%
  dplyr::mutate(Matrix = factor(Matrix, levels = Samples_levels)) %>%
  dplyr::arrange(Matrix)

for (i in seq_along(recovery_range$Matrix)) {
  recovery_range$Matrix_label[i] <- paste(recovery_range$`Euclidean Distance`[i], " (", recovery_range$mean_lab[i], ") ", recovery_range$Matrix[i], sep = "")
  recovery_range$mean_lab_let[i] <- paste("(", letters[i], ") ", recovery_range$mean_lab[i], sep = "")
  recovery_range$dist.recov[i] <- paste(recovery_range$`Euclidean Distance`[i], " (", recovery_range$mean_lab[i], ")", sep = "")
}

AP_cluster_plot$recovery_label <- NA
AP_cluster_plot$matrix_label <- NA
AP_cluster_plot$dist.recov <- NA

for (i in seq_along(AP_cluster_plot$Sample.Name)) {
  if(AP_cluster_plot$Sample.Name[i] %in% curve_points) {
    AP_cluster_plot$recovery_label[i] <- ""
    AP_cluster_plot$dist.recov[i] <- ""
    AP_cluster_plot$matrix_label[i] <- as.character(AP_cluster_plot$Sample.Name[i])
  }
  
  else {
    AP_cluster_plot$recovery_label[i] <- recovery_range$mean_lab_let[as.character(recovery_range$Matrix) == as.character(AP_cluster_plot$Sample.Name[i])]
    AP_cluster_plot$dist.recov[i] <- recovery_range$dist.recov[as.character(recovery_range$Matrix) == as.character(AP_cluster_plot$Sample.Name[i])]
    AP_cluster_plot$matrix_label[i] <- recovery_range$Matrix_label[as.character(recovery_range$Matrix) == as.character(AP_cluster_plot$Sample.Name[i])]
  }
}

samples_levels_let <- NULL

for (i in seq_along(Samples_levels)) {
  samples_levels_let[i] <- paste("(", letters[i],"). ", Samples[i], sep = "")
}

AP_cluster_plot$Sample.Name <- as.character(AP_cluster_plot$Sample.Name)

for (i in seq_along(AP_cluster_plot$Sample.Name)) {
  if(AP_cluster_plot$Sample.Name[i] %in% curve_points) {
    AP_cluster_plot$Sample.Name[i] <- as.character(AP_cluster_plot$Sample.Name[i])
  }
  
  else {
    AP_cluster_plot$Sample.Name[i] <- paste("(",letters[i], "). ", as.character(AP_cluster_plot$Sample.Name[i]), sep = "")
    AP_cluster_plot$dist.recov[i] <- paste("(",letters[i], "). ", AP_cluster_plot$dist.recov[i], sep = "")
  }
}

AP_cluster_plot_filt <- AP_cluster_plot %>%
  dplyr::mutate(matrix_factor = case_when(
    Sample.Name %in% curve_points ~ "Calibration curve points",
    !Sample.Name %in% curve_points ~ Sample.Name
  )) %>%
  dplyr::mutate(matrix_factor = factor(matrix_factor, levels = c("Calibration curve points", samples_levels_let)))

# Figure 6 ====

AP_PCs_with_distance <- AP_cluster_plot_filt %>%
  ggplot(aes(PC1, PC2, label = dist.recov)) +
  geom_point(aes(fill = matrix_factor,
                 shape = matrix_factor,
                 size = matrix_factor),
             color = "black",
             alpha = 0.75) +
  scale_shape_manual(values = c(4,
                                21, 21, 21, 21,
                                22, 22,
                                25, 25,
                                24, 24, 24, 24,
                                23, 23),
                     labels = c("Calibration curve points", samples_levels_let)) +
  scale_size_manual(values = c(1,
                               3, 3, 3, 3,
                               3, 3,
                               3, 3,
                               3, 3, 3, 3,
                               3, 3),
                    labels = c("Calibration curve points", samples_levels_let)) +
  scale_fill_manual(values = c("black",
                               viridis(4, alpha = 0.9)[4], viridis(4, alpha = 0.9)[3], viridis(4, alpha = 0.9)[2], viridis(4, alpha = 0.9)[1],
                               viridis(4, alpha = 0.9)[4], viridis(4, alpha = 0.9)[1],
                               viridis(4, alpha = 0.9)[4], viridis(4, alpha = 0.9)[1],
                               viridis(4, alpha = 0.9)[4], viridis(4, alpha = 0.9)[3], viridis(4, alpha = 0.9)[2], viridis(4, alpha = 0.9)[1],
                               viridis(4, alpha = 0.9)[4], viridis(4, alpha = 0.9)[1]),
                    labels = c("Calibration curve points", samples_levels_let)) +
  geom_segment(aes(x = PC1, y = PC2, xend = Exemp_x, yend = Exemp_y,
                   group = Group),
               alpha = 0.5) +
  geom_text_repel(segment.alpha = 0.5,
                  size = 2.5) +
  guides(size = "none") +
  theme_bw() +
  theme(text = element_text(family = "Times",
                            face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 14),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.key.size = unit(0.15, "inches"),
        panel.grid = element_blank()) +
  labs(x = "PC 1",
       y = "PC 2",
       fill = "Matrix \n      ",
       shape = "Matrix \n      ")
AP_PCs_with_distance

ggsave(paste(getwd(), "/manuscript/color_theme/Fig.6.tiff", sep = ""), plot = AP_PCs_with_distance, units = "in", width = 6, height = 4, dpi = 1000, compression = "lzw")
dev.off()
