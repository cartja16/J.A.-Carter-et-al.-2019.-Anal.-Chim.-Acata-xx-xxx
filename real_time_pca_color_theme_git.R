# Notes ===

# The following provides the code for PCA and AP Clustering to produce Figures 1, SI1-SI3

# Clear worksapce
rm(list=ls()) 
dev.off()

# load pacakges
library(tidyverse)
library(broom)
library(viridis)
library(ggrepel)
library(apcluster)
library(RColorBrewer)
library(extrafont)
library(readxl)
library(gridExtra)

windowsFonts(Times=windowsFont("TT Times New Roman"))
par(family = "Times")

analytes <- c("Cd.214.439", "Co.238.892", "Cr.267.716", "Pb.220.353")
intstd <- c("Y.360.074", "Y.371.029", "Y.377.433")
curve_points <- c("0 ppb", "100 ppb", "200 ppb", "500 ppb", "1000 ppb")

# Read in and format data ====

matrix <- read.csv("real_time_Na_Carb_matrix_091118_fit.csv") # here should be the input for the csv, change for the matrix in question
matrix_curve <- read.csv("curve_LOD_091118_Fit.csv") %>%
  dplyr::rename(`Sample Name` = Label) %>%
  dplyr::mutate_at(vars(Element), funs(gsub(" ", ".", ., fixed = TRUE))) %>%
  dplyr::mutate(`Sample Name` = case_when(grepl("Blank", `Sample Name`) ~ "0 ppb",
                                          grepl("Standard 1", `Sample Name`) ~ "100 ppb",
                                          grepl("Standard 2", `Sample Name`) ~ "200 ppb",
                                          grepl("Standard 3", `Sample Name`) ~ "500 ppb",
                                          grepl("Standard 4", `Sample Name`) ~ "1000 ppb")
  ) %>%
  tidyr::gather(Replicate, Intensity, -`Sample Name`, -Element) %>%
  tidyr::spread(Element, Intensity) %>%
  dplyr::group_by(`Sample Name`) %>%
  dplyr::summarize_at(vars(-Replicate), funs(mean(., na.rm=TRUE))) %>%
  dplyr::ungroup() 
  
BGVars_intstd <- matrix %>%
  dplyr::select(-Time.min, -analytes) %>%
  colnames()

BGVars <- matrix %>%
  dplyr::select(-Time.min, -analytes, -intstd) %>%
  colnames()

time_points <- matrix %>%
  dplyr::select(Time.min) %>%
  dplyr::mutate(Time.min = as.character(Time.min)) %>%
  dplyr::pull()

# PCA stats whole ====

PCA_curve_mean <- function(prcomp_object, OES_norm_wide, curve_character_string, sample_character_string) {
  tmp_avg <- as.tibble(prcomp_object$x) %>%
    dplyr::filter(OES_norm_wide$`Sample Name` %in% curve_character_string) %>%
    dplyr::summarize_at(vars(PC1, PC2), funs(mean(., na.rm=TRUE))) %>%
    as.numeric()
  
  tmp_samps <- as.tibble(prcomp_object$x) %>%
    dplyr::filter(OES_norm_wide$`Sample Name` %in% sample_character_string) %>%
    dplyr::select(PC1, PC2) 
  
  tmp_distances <- lst()
  for (i in seq_along(sample_character_string)) {
    tmp_distances[i] <- sqrt((tmp_samps[i,1] - tmp_avg[1])^2 + (tmp_samps[i,2] - tmp_avg[2])^2)
  }
  names(tmp_distances) <- sample_character_string
  
  distances_out <- tibble(seq_along(tmp_distances))
  for (j in seq_along(tmp_distances)) {
    distances_out[j,1] <- tmp_distances[[j]]
    colnames(distances_out) <- "Euclidean_distance"
  }
  
  distances_out <- distances_out %>%
    dplyr::mutate(`Sample Name` = sample_character_string) %>%
    dplyr::select(`Sample Name`, Euclidean_distance)
  return(distances_out)
}

matrix_BG_norm <- matrix %>%
  dplyr::rename(`Sample Name` = Time.min) %>%
  dplyr::mutate(`Sample Name` = as.character(`Sample Name`)) %>%
  dplyr::bind_rows(matrix_curve) %>%
  dplyr::select(-analytes, -intstd) %>%
  tidyr::gather(Element, Intensity, -`Sample Name`) %>%
  dplyr::group_by(Element) %>%
  dplyr::mutate(normalized = Intensity/max(Intensity[`Sample Name` %in% curve_points])) %>%
  dplyr::select(-Intensity) %>%
  tidyr::spread(Element, normalized) %>%
  dplyr::mutate(`Sample Name` = as.character(`Sample Name`))

matrix_pca <- matrix_BG_norm %>%
  dplyr::select(-`Sample Name`) %>%
  scale(center = TRUE, scale = TRUE) %>%
  prcomp()

matrix_dist_Euclidean_calc <- PCA_curve_mean(matrix_pca, matrix_BG_norm, curve_points, time_points)

matrix_dist_tmp <- matrix_BG_norm %>%
  filter(`Sample Name` %in% curve_points) %>%
  dplyr::mutate(Euclidean_distance = "") 

matrix_dist <- matrix_BG_norm %>%
  dplyr::filter(`Sample Name` %in% time_points) %>%
  dplyr::mutate(Euclidean_distance = case_when(
    `Sample Name` == matrix_dist_Euclidean_calc$`Sample Name` ~ as.character(signif(matrix_dist_Euclidean_calc$Euclidean_distance, digits = 2))
  )) %>%
  dplyr::bind_rows(matrix_dist_tmp) %>%
  dplyr::mutate(`Sample Name` = factor(`Sample Name`, levels = c(curve_points, time_points),ordered = TRUE))

# For correct row order

matrix_dist_pca <- matrix_dist %>%
  dplyr::select(-`Sample Name`, -Euclidean_distance) %>%
  scale(center = TRUE, scale = TRUE) %>%
  prcomp()

matrix_dist <- matrix_dist %>%
  dplyr::bind_cols(as.tibble(matrix_dist_pca$x))

# Figure 3

Euclid_bar_whole <- matrix_dist %>%
  dplyr::filter(!`Sample Name` %in% curve_points) %>%
  dplyr::select(`Sample Name`, Euclidean_distance) %>%
  dplyr::mutate(`Sample Name` = as.numeric(time_points)) %>%
  dplyr::arrange(`Sample Name`) %>%
  dplyr::mutate(Euclidean_distance = as.numeric(Euclidean_distance)) %>%
  ggplot(aes(`Sample Name`, Euclidean_distance)) +
  geom_bar(stat = "identity",
           color = "blue",
           fill = "white") +
  xlab("Time (min)") +
  ylab("Euclidean Distance") +
  geom_smooth(se = FALSE,
              color = "red") +
  theme_bw() +
  theme(text = element_text(family = "Times",
                            face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        panel.grid = element_blank())
Euclid_bar_whole
#ggsave(paste(getwd(), "/manuscript/color_theme/Fig.3.tiff", sep =""), plot = Euclid_bar_whole, units = "in", width = 6, height = 4, compression = "lzw")
#dev.off()

# AP Cluster whole ====

AP_cluster_whole <- matrix_dist %>%
  dplyr::select(-BGVars, -`Sample Name`, -Euclidean_distance) 

AP_cluster_sim_whole <- negDistMat(AP_cluster_whole, r=2)
AP_cluster_model_whole <- apcluster(AP_cluster_sim_whole, details = TRUE)

AP_cluster_plot_whole <- matrix_dist

Exemplars_ind_whole <- AP_cluster_model_whole@exemplars

AP_cluster_plot_whole$Exemplars <- "Non_Exemplar"

AP_cluster_plot_whole$Exemplars[Exemplars_ind_whole] <- "Exemplar"

AP_cluster_plot_whole <- AP_cluster_plot_whole %>%
  dplyr::mutate(Group = factor(AP_cluster_model_whole@idx, levels = unique(AP_cluster_model_whole@idx), ordered = TRUE)) %>%
  dplyr::mutate(Exemplars = factor(Exemplars, levels = c("Non_Exemplar", "Exemplar"))) %>%
  dplyr::mutate(`Sample Name` = factor(`Sample Name`, levels = c(curve_points, time_points, ordered = TRUE))) %>%
  dplyr::group_by(Group) %>%
  dplyr::mutate(Exemp_x = PC1[Exemplars == "Exemplar"],
                Exemp_y = PC2[Exemplars == "Exemplar"]) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(sample_chr = as.character(`Sample Name`)) %>%
  dplyr::mutate(color_cont = case_when(
    sample_chr %in% curve_points ~ 0,
    !sample_chr %in% curve_points ~ as.numeric(sample_chr)
  )) %>%
  dplyr::mutate(Euclidean_distance_num = case_when(
    sample_chr %in% curve_points ~ 0,
    !sample_chr %in% curve_points ~ as.numeric(Euclidean_distance)
  ))

max_whole <- max(AP_cluster_plot_whole$Euclidean_distance_num)

# Figures 2, SI4-SI6 ====

AP_PCs_whole <- ggplot(AP_cluster_plot_whole, aes(PC1, PC2)) +
  geom_point(aes(fill = Euclidean_distance_num),
             color = "black",
             stroke = 1,
             size = 4.5,
             shape = 21,
             alpha = 0.75) +
  scale_fill_viridis_c(guide = guide_colorbar(frame.colour = "black",
                                             frame.linewidth = 1,
                                             ticks.colour = "white",
                                             ticks.linewidth = 1),
                       direction = -1,
                       limits = c(0, plyr::round_any(max_whole, 0.5, ceiling)),
                       breaks = c(round((plyr::round_any(max_whole, 0.5, ceiling)/4)*1, digits = 2),
                                  round((plyr::round_any(max_whole, 0.5, ceiling)/4)*2, digits = 2),
                                  round((plyr::round_any(max_whole, 0.5, ceiling)/4)*3, digits = 2),
                                  round((plyr::round_any(max_whole, 0.5, ceiling)/4)*4, digits = 2))) +
  geom_segment(aes(x = PC1, y = PC2, xend = Exemp_x, yend = Exemp_y,
                   group = Group),
               alpha = 0.5) +
  labs(fill = "Euclidean distances from \n calibration curve points \n              ") +
  theme_bw() +
  theme(text = element_text(family = "Times",
                            face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.key.size = unit(0.15, "inches"),
        panel.grid = element_blank())
AP_PCs_whole
#ggsave(paste(getwd(), "/manuscript/color_theme/", "Fig.2.tiff", sep =""), plot = AP_PCs_whole, units = "in", width = 6, height = 4, compression = "lzw")

# PCA stats iterative ====

pca_build <- lst()
pca_real_time <- tibble()
Euclid_tib <- tibble()

curve_tmp <- matrix_BG_norm %>%
  dplyr::filter(`Sample Name` %in% curve_points)

for (i in seq_along(time_points)) {
  tib_tmp <- matrix_BG_norm %>%
    dplyr::filter(`Sample Name` %in% time_points) %>%
    dplyr::slice(i) %>%
    dplyr::bind_rows(curve_tmp)
  
  pca_tmp <- tib_tmp %>%
    dplyr::select(-`Sample Name`) %>%
    scale(center = TRUE, scale = TRUE) %>%
    prcomp()
  
  Euclid_tib <- rbind(Euclid_tib, PCA_curve_mean(pca_tmp, tib_tmp, curve_points, time_points[i]))
  
  pca_build[[i]] <- as.tibble(pca_tmp$x) %>%
    dplyr::bind_cols(tib_tmp)
  
  pca_real_time <- rbind(pca_real_time, pca_build[[i]])
}

pca_real_time_curve <- pca_real_time %>%
  dplyr::filter(`Sample Name` %in% curve_points) %>%
  dplyr::mutate(Euclidean_distance = "")

pca_real_time <- pca_real_time %>%
  dplyr::filter(`Sample Name` %in% time_points) %>%
  dplyr::mutate(Euclidean_distance = as.character(signif(Euclid_tib$Euclidean_distance, digits = 3))) %>%
  dplyr::bind_rows(pca_real_time_curve)

# AP Cluster iterative ====

AP_cluster <- pca_real_time %>%
  dplyr::select(-BGVars, -`Sample Name`, -Euclidean_distance) 

AP_cluster_sim <- negDistMat(AP_cluster, r=2)
AP_cluster_model <- apcluster(AP_cluster_sim, details = TRUE)

AP_cluster_plot <- pca_real_time

Exemplars_ind <- AP_cluster_model@exemplars

AP_cluster_plot$Exemplars <- "Non_Exemplar"

AP_cluster_plot$Exemplars[Exemplars_ind] <- "Exemplar"

AP_cluster_plot <- AP_cluster_plot %>%
  dplyr::mutate(Group = factor(AP_cluster_model@idx, levels = unique(AP_cluster_model@idx), ordered = TRUE)) %>%
  dplyr::mutate(Exemplars = factor(Exemplars, levels = c("Non_Exemplar", "Exemplar"))) %>%
  dplyr::mutate(`Sample Name` = factor(`Sample Name`, levels = c(curve_points, time_points), ordered = TRUE)) %>%
  dplyr::group_by(Group) %>%
  dplyr::mutate(Exemp_x = PC1[Exemplars == "Exemplar"],
                Exemp_y = PC2[Exemplars == "Exemplar"]) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(sample_chr = as.character(`Sample Name`)) %>%
  dplyr::mutate(color_cont = case_when(
    sample_chr %in% curve_points ~ 0,
    !sample_chr %in% curve_points ~ as.numeric(sample_chr)
  )) %>%
  dplyr::mutate(Euclidean_distance_num = case_when(
    sample_chr %in% curve_points ~ 0,
    !sample_chr %in% curve_points ~ as.numeric(Euclidean_distance)
  ))

max_iter <- max(AP_cluster_plot$Euclidean_distance_num)

round(plyr::round_any(max_iter, 0.5, ceiling)/1, digits = 2)

# Figures SI7 - SI11 ====

AP_PCs <- ggplot(AP_cluster_plot, aes(PC1, PC2)) +
  geom_point(aes(fill = Euclidean_distance_num),
             color = "black",
             stroke = 1,
             size = 4.5,
             shape = 21,
             alpha = 0.75) +
  scale_fill_viridis_c(guide = guide_colorbar(frame.colour = "black",
                                              frame.linewidth = 1,
                                              ticks.colour = "white",
                                              ticks.linewidth = 1),
                       direction = -1,
                       limits = c(0, ceiling(max_iter)),
                       breaks = c(round((plyr::round_any(max_iter, 0.5, ceiling)/4)*1, digits = 2),
                                  round((plyr::round_any(max_iter, 0.5, ceiling)/4)*2, digits = 2),
                                  round((plyr::round_any(max_iter, 0.5, ceiling)/4)*3, digits = 2),
                                  round((plyr::round_any(max_iter, 0.5, ceiling)/4)*4, digits = 2))) +
  geom_segment(aes(x = PC1, y = PC2, xend = Exemp_x, yend = Exemp_y,
                   group = Group),
               alpha = 0.5) +
  labs(fill = "Euclidean distances from \n calibration curve points \n              ") +
  theme_bw() +
  theme(text = element_text(family = "Times",
                            face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.key.size = unit(0.15, "inches"),
        panel.grid = element_blank())
AP_PCs

#ggsave(paste(getwd(), "/manuscript/color_theme/", "Fig.SI.7.tiff", sep =""), plot = AP_PCs, units = "in", width = 6, height = 4, compression = "lzw")
#dev.off()
