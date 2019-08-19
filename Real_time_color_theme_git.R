# Notes ====

# The following provides the code for Figures 1, SI1-SI3

# Clear workspace
rm(list=ls())
dev.off()

# load packages
library(readxl)
library(gridExtra)
library(extrafont)
library(tidyverse)
library(broom)
library(viridis)
library(RColorBrewer)

windowsFonts(Times=windowsFont("TT Times New Roman"))
par(family = "Times")

# Read in and format data ====

analytes <- c("Cd.214.439", "Co.238.892", "Cr.267.716", "Pb.220.353")
intstd <- c("Y.360.074", "Y.371.029", "Y.377.433")
analytes_Cd <- c("Co.238.892", "Cr.267.716", "Pb.220.353")
analytes_Co <- c("Cd.214.439", "Cr.267.716", "Pb.220.353")
analytes_Cr <- c("Cd.214.439", "Co.238.892", "Pb.220.353")
analytes_Pb <- c("Cd.214.439", "Co.238.892", "Cr.267.716")

matrix <- read.csv("real_time_Na_Carb_matrix_091118_fit.csv")
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

# Functions to process and format data for plotting ====

BG_ratio_product <- function(analyte_vector_rawIntensities, BG_string, df_with_all_intensities, analytes_selection_string) {
  correction_out <- list()
  for (i in seq_along(BG_string)) {
    correction_out[[i]] <- analyte_vector_rawIntensities/df_with_all_intensities[BG_string[i]]
    names(correction_out)[i] <- BG_string[i]
    names(correction_out)[i] <- paste("ratio", names(correction_out)[i], sep = ".")
  }
  
  l <- length(BG_string)
  
  for (j in seq_along(BG_string)) {
    correction_out[[j + l]] <- analyte_vector_rawIntensities*df_with_all_intensities[BG_string[j]]
    names(correction_out)[j + l] <- BG_string[j]
    names(correction_out)[j + l] <- paste("product", names(correction_out)[j + l], sep = ".")
  }
  
  correction_named <- as.data.frame(correction_out)
  colnames(correction_named) <- names(correction_out)
  
  correction_named <- cbind(df_with_all_intensities, correction_named)
  correction_named <- as_tibble(correction_named) %>%
    dplyr::select(-analytes_selection_string)
  return(correction_named)
}

samples_calculation <- function(samples_intensity_df, analyte_stats_df, analyte_character, curve_tibble) {
  curve_tmp <- curve_tibble %>%
    dplyr::filter(conc == 0)
  
  for (i in seq_along(samples_intensity_df$conc)) {
    samples_intensity_df$conc[i] = ((samples_intensity_df$Intensity[i] - curve_tmp$Intensity[curve_tmp$Element == samples_intensity_df$Element[i]]) - analyte_stats_df[which(analyte_stats_df$Element == samples_intensity_df$Element[i]),]["estimate"][1,]) / analyte_stats_df[which(analyte_stats_df$Element == samples_intensity_df$Element[i]),]["estimate"][2,]
  }
  samples_intensity_df$conc <- as.numeric(samples_intensity_df$conc)
  concs_tmp <- samples_intensity_df %>%
    dplyr::select(-Intensity) %>%
    tidyr::spread(Element, conc)
  
  tmp_stats <- concs_tmp %>%
    dplyr::summarize(mean_raw = mean(!!analyte_character),
                     max_raw = max(!!analyte_character),
                     min_raw = min(!!analyte_character),
                     sd_raw = sd(!!analyte_character),
                     rsd_raw = (sd(!!analyte_character)/mean(!!analyte_character)) * 100
                     )
  
  concs_out <- concs_tmp %>%
    dplyr::mutate(
      mean_raw = tmp_stats$mean_raw,
      max_raw = tmp_stats$max_raw,
      min_raw = tmp_stats$min_raw,
      sd_raw = tmp_stats$sd_raw,
      rsd_raw = tmp_stats$rsd_raw
      )
  
  return(concs_out)
}

curve_build <- function(analyte_vector_intensities, BG_vars_character, curve_tibble, analytes_vector, net = TRUE) {
  if(net) {
    curve_out <- BG_ratio_product(analyte_vector_intensities, BG_vars_character, curve_tibble, analytes_vector) %>%
      dplyr::select(-BGVars_intstd) %>%
      dplyr::mutate(conc = case_when(
        `Sample Name` == "0 ppb" ~ 0,
        `Sample Name` == "100 ppb" ~ 100,
        `Sample Name` == "200 ppb" ~ 200,
        `Sample Name` == "500 ppb" ~ 500,
        `Sample Name` == "1000 ppb" ~ 1000
      )) %>%
      tidyr::gather(Element, Intensity, -`Sample Name`, -conc) %>%
      dplyr::group_by(Element) %>%
      dplyr::mutate(
        Intensity = Intensity - Intensity[`Sample Name` == "0 ppb"]
      ) %>%
      dplyr::ungroup()
    return(curve_out)
  }
  if(!net) {
    curve_out <- BG_ratio_product(analyte_vector_intensities, BG_vars_character, curve_tibble, analytes_vector) %>%
      dplyr::select(-BGVars_intstd) %>%
      dplyr::mutate(conc = case_when(
        `Sample Name` == "0 ppb" ~ 0,
        `Sample Name` == "100 ppb" ~ 100,
        `Sample Name` == "200 ppb" ~ 200,
        `Sample Name` == "500 ppb" ~ 500,
        `Sample Name` == "1000 ppb" ~ 1000
      )) %>%
      tidyr::gather(Element, Intensity, -`Sample Name`, -conc)
    return(curve_out)
  }
}

sample_build <- function(analyte_vector_intensities, BG_vars_character, sample_tibble, analytes_vector) {
  sample_correction_out <- BG_ratio_product(analyte_vector_intensities, BG_vars_character, sample_tibble, analytes_vector) %>%
    dplyr::mutate(conc = 0) %>%
    dplyr::select(-BGVars_intstd) %>%
    tidyr::gather(Element, Intensity, -Time.min, -conc)
  return(sample_correction_out)
}

analyte_curve_stats <- function(curve_df) {
  df_out <- curve_df %>%
    tidyr::nest(-Element, `Sample Name`) %>% 
    dplyr::mutate(
      fit = purrr::map(data, ~ lm(Intensity ~ conc, data = .x)),
      tidied = purrr::map(fit, tidy),
      glanced = purrr::map(fit, glance)
    ) %>% 
    unnest(tidied)
  return(df_out)
}

time_plot_gg <- function(BG_correction_tibble, signal_character, title) {
  
  g <- ggplot(BG_correction_tibble, aes(Time.min, !!signal_character)) +
    geom_point(fill = brewer.pal(3, "Greys")[1],
               color = "blue",
               size = 4.5,
               stroke = 1,
               shape = 21,
               alpha = 0.75) +
    geom_smooth(se = FALSE,
                color = "red",
                size = 1,
                alpha = 0.75,
                method = "loess") + 
    geom_hline(yintercept = c(500, 450, 550, 400, 600),
               linetype = c("dashed", "solid", "solid", "longdash", "longdash")) +
    scale_y_continuous(limits = c(150,700),
                       breaks=seq(150,700,50)) +
    ggtitle(title) +
    xlab("Time (minutes)") +
    theme_bw() +
    theme(text = element_text(face = "bold",
                              family = "Times"),
          axis.text = element_text(size = 14),
          title = element_text(size = 14))
  return(g)
}

# Data processing ====

matrix_BG_correction_Cd <- sample_build(matrix$Cd.214.439, BGVars_intstd, matrix, analytes_Cd) 
matrix_BG_correction_Co <- sample_build(matrix$Co.238.892, BGVars_intstd, matrix, analytes_Co) 
matrix_BG_correction_Cr <- sample_build(matrix$Cr.267.716, BGVars_intstd, matrix, analytes_Cr)
matrix_BG_correction_Pb <- sample_build(matrix$Pb.220.353, BGVars_intstd, matrix, analytes_Pb) 

matrix_BG_curve_Cd_net <- curve_build(matrix_curve$Cd.214.439, BGVars_intstd, matrix_curve, analytes_Cd, net = TRUE)
matrix_BG_curve_Cd <- curve_build(matrix_curve$Cd.214.439, BGVars_intstd, matrix_curve, analytes_Cd, net = FALSE)

matrix_BG_curve_Co_net <- curve_build(matrix_curve$Co.238.892, BGVars_intstd, matrix_curve, analytes_Co, net = TRUE)
matrix_BG_curve_Co <- curve_build(matrix_curve$Co.238.892, BGVars_intstd, matrix_curve, analytes_Co, net = FALSE)

matrix_BG_curve_Cr_net <- curve_build(matrix_curve$Cr.267.716, BGVars_intstd, matrix_curve, analytes_Cr, net = TRUE)
matrix_BG_curve_Cr <- curve_build(matrix_curve$Cr.267.716, BGVars_intstd, matrix_curve, analytes_Cr, net = FALSE)

matrix_BG_curve_Pb_net <- curve_build(matrix_curve$Pb.220.353, BGVars_intstd, matrix_curve, analytes_Pb, net = TRUE)
matrix_BG_curve_Pb <- curve_build(matrix_curve$Pb.220.353, BGVars_intstd, matrix_curve, analytes_Pb, net = FALSE)

matrix_curve_stats_Cd <- analyte_curve_stats(matrix_BG_curve_Cd_net)
matrix_curve_stats_Co <- analyte_curve_stats(matrix_BG_curve_Co_net)
matrix_curve_stats_Cr <- analyte_curve_stats(matrix_BG_curve_Cr_net)
matrix_curve_stats_Pb <- analyte_curve_stats(matrix_BG_curve_Pb_net)

matrix_Cd_concs <- samples_calculation(matrix_BG_correction_Cd, matrix_curve_stats_Cd, quo(Cd.214.439), matrix_BG_curve_Cd)
matrix_Co_concs <- samples_calculation(matrix_BG_correction_Co, matrix_curve_stats_Co, quo(Co.238.892), matrix_BG_curve_Co)
matrix_Cr_concs <- samples_calculation(matrix_BG_correction_Cr, matrix_curve_stats_Cr, quo(Cr.267.716), matrix_BG_curve_Cr)
matrix_Pb_concs <- samples_calculation(matrix_BG_correction_Pb, matrix_curve_stats_Pb, quo(Pb.220.353), matrix_BG_curve_Pb)

Cd_raw <- time_plot_gg(matrix_Cd_concs, quo(Cd.214.439), paste("(A) RSD =", as.character(signif(matrix_Cd_concs$rsd_raw, digits = 2)), "%"))

Co_raw <- time_plot_gg(matrix_Co_concs, quo(Co.238.892), paste("(B) RSD =", as.character(signif(matrix_Co_concs$rsd_raw, digits = 2)), "%"))

Cr_raw <- time_plot_gg(matrix_Cr_concs, quo(Cr.267.716), paste("(C) RSD =", as.character(signif(matrix_Cr_concs$rsd_raw, digits = 2)), "%"))

Pb_raw <- time_plot_gg(matrix_Pb_concs, quo(Pb.220.353), paste("(D) RSD =", as.character(signif(matrix_Pb_concs$rsd_raw, digits = 2)), "%"))

# Figures 1, SI1-SI3 ====

raw_gg <- arrangeGrob(grobs = list(Cd_raw, Co_raw, Cr_raw, Pb_raw),
                      nrow = 2, ncol = 2)

ggsave(paste(getwd(), "/manuscript/color_theme/test.tiff", sep = ""), plot = raw_gg, units = "in", width = 10, height = 6, dpi = 1000, compression = "lzw")
dev.off()
