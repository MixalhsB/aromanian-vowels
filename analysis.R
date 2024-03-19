library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(emuR)


# Some functions:

get_speaker_specific_MD <- function(dataset, speaker) {
  current_subset <- subset(dataset, Speaker == speaker)
  MD <- mahalanobis(current_subset[, c(11,12,13)], # c(11, 12, 13) = F1, F2, F3
                    colMeans(current_subset[, c(11,12,13)]),
                    cov(current_subset[, c(11,12,13)]))
  return(MD)
}

get_speaker_vowel_specific_MD <- function(dataset, speaker, vowel) {
  current_subset <- subset(dataset, Speaker == speaker & VowelType == vowel)
  MD <- mahalanobis(current_subset[, c(11,12,13)], # c(11, 12, 13) = F1, F2, F3
                    colMeans(current_subset[, c(11,12,13)]),
                    cov(current_subset[, c(11,12,13)]))
  return(MD)
}


# Initial preparation of datasets:

aromanian <- read.csv('./aromanian_targets.csv')
all_vowels <- read.csv('./aromanian_all_raw.csv')

aromanian$Sex <- ifelse(substring(aromanian$Speaker, 1, 1) == "W", "Female", "Male")
all_vowels$Sex <- ifelse(substring(all_vowels$Speaker, 1, 1) == "W", "Female", "Male")

aromanian$InitialID <- sapply(sapply(seq(nrow(aromanian)), toString), function(x) paste0("T", str_pad(x, 4, pad = "0")))
all_vowels$InitialID <- sapply(sapply(seq(nrow(all_vowels)), toString), function(x) paste0("V", str_pad(x, 5, pad = "0")))

aromanian <- mutate(aromanian, VowelType = as.factor(VowelType),
                    Dialect = as.factor(Dialect), Speaker = as.factor(Speaker),
                    POS = as.factor(POS), Stress = as.factor(Stress))
all_vowels <- mutate(all_vowels, VowelType = as.factor(substring(str_extract(VowelInWord, "\\(.+\\)"), 2, 2)),
                     Dialect = as.factor(Dialect), Speaker = as.factor(Speaker))

aromanian <- mutate(aromanian, F1 = bark(F1), F2 = bark(F2), F3 = bark(F3))
all_vowels <- mutate(all_vowels, F1 = bark(F1), F2 = bark(F2), F3 = bark(F3))


# MD outlier exclusion:

already_filtered <- FALSE
plots <- list()
other_i <- 0

if (!already_filtered) {
  
  ## step 1, for targets:
  aromanian <- aromanian[order(aromanian$VowelType), ] # now sorted by vowel (and speaker)
  MD_moderate_thresholds <- c()
  MD_extreme_thresholds <- c()
  complete_MD_vector <- c()
  for (vowel in unique(aromanian$VowelType)) {
    for (speaker in unique(aromanian$Speaker)) {
      MD <- get_speaker_vowel_specific_MD(aromanian, speaker, vowel)
      MD_moderate_thresholds <- c(MD_moderate_thresholds, 7.82)
      MD_extreme_thresholds <- c(MD_extreme_thresholds, 7.82)
      complete_MD_vector <- c(complete_MD_vector, MD)
    }
  }
  aromanian$MDBySpeaker <- complete_MD_vector
  
  ## step 1, for all:
  all_vowels <- all_vowels[order(all_vowels$VowelType), ] # now sorted by vowel (and speaker)
  MD_moderate_thresholds <- c()
  MD_extreme_thresholds <- c()
  complete_MD_vector <- c()
  for (vowel in unique(all_vowels$VowelType)) {
    for (speaker in unique(all_vowels$Speaker)) {
      MD <- get_speaker_vowel_specific_MD(all_vowels, speaker, vowel)
      MD_moderate_thresholds <- c(MD_moderate_thresholds, 7.82)
      MD_extreme_thresholds <- c(MD_extreme_thresholds, 7.82)
      complete_MD_vector <- c(complete_MD_vector, MD)
    }
  }
  all_vowels$MDBySpeaker <- complete_MD_vector
  
  # step 2, for targets:
  complete_IsModOut_vector <- c()
  complete_IsExtOut_vector <- c()
  k <- 0
  for (i in 1:length(unique(aromanian$VowelType))) {
    for (j in 1:length(unique(aromanian$Speaker))) {
      k <- k + 1
      vowel <- unique(aromanian$VowelType)[i]
      speaker <- unique(aromanian$Speaker)[j]
      moderate_threshold <- MD_moderate_thresholds[k]
      extreme_threshold <- MD_extreme_thresholds[k]
      
      current_subset <- subset(aromanian, Speaker == speaker & VowelType == vowel)
      
      complete_IsModOut_vector <- c(complete_IsModOut_vector, current_subset$MDBySpeaker > moderate_threshold) 
      
      complete_IsExtOut_vector <- c(complete_IsExtOut_vector, current_subset$MDBySpeaker > extreme_threshold)
    }
  }
  aromanian <- mutate(aromanian, IsModerateOut = complete_IsModOut_vector)
  aromanian <- mutate(aromanian, IsExtremeOut = complete_IsExtOut_vector)
  aromanian <- aromanian[order(aromanian$InitialID), ] # revert sorting
  
  # step 2, for all:
  complete_IsModOut_vector <- c()
  complete_IsExtOut_vector <- c()
  k <- 0
  for (i in 1:length(unique(all_vowels$VowelType))) {
    for (j in 1:length(unique(all_vowels$Speaker))) {
      k <- k + 1
      vowel <- unique(all_vowels$VowelType)[i]
      speaker <- unique(all_vowels$Speaker)[j]
      moderate_threshold <- MD_moderate_thresholds[k]
      extreme_threshold <- MD_extreme_thresholds[k]
      
      current_subset <- subset(all_vowels, Speaker == speaker & VowelType == vowel)
      
      complete_IsModOut_vector <- c(complete_IsModOut_vector, current_subset$MDBySpeaker > moderate_threshold)
      
      complete_IsExtOut_vector <- c(complete_IsExtOut_vector, current_subset$MDBySpeaker > extreme_threshold)
    }
  }
  all_vowels <- mutate(all_vowels, IsModerateOut = complete_IsModOut_vector)
  all_vowels <- mutate(all_vowels, IsExtremeOut = complete_IsExtOut_vector)
  all_vowels <- all_vowels[order(all_vowels$InitialID), ] # revert sorting
  
  # Filtering out for both:
  
  aromanian <- subset(aromanian, IsExtremeOut == FALSE)
  aromanian <- subset(aromanian, select = -c(IsExtremeOut, IsModerateOut))
  
  all_vowels <- subset(all_vowels, IsExtremeOut == FALSE)
  all_vowels <- subset(all_vowels, select = -c(IsExtremeOut, IsModerateOut))
  
  already_filtered = TRUE
}


### Computing Pillai scores

factors <- c("VowelType")

complete_cmb_list <- list()
for (speaker in unique(aromanian$Speaker)) {
  current_subset <- subset(aromanian, Speaker == speaker)
  cmb <- list()

  cmb[[factors[1]]] <- manova(cbind(F1, F2, F3) ~ VowelType, data = current_subset)

  complete_cmb_list[[speaker]] <- cmb
}

target_pillai_results <- c()
for (speaker in unique(aromanian$Speaker)) {
  for (fct in factors) {
    summary_mnv <- summary(complete_cmb_list[[speaker]][[fct]])
    pillai <- summary_mnv$stats[, "Pillai"][fct]
    target_pillai_results <- c(target_pillai_results, pillai)
  }
}

target_pillai_matrix <- matrix(target_pillai_results, nrow = length(factors), ncol = length(unique(aromanian$Speaker)))

colnames(target_pillai_matrix) <- unique(aromanian$Speaker)

target_pillai_rownames <- c()
for (fct in factors) {
  target_pillai_rownames <- c(target_pillai_rownames, paste(fct, "Pillai"))
}
rownames(target_pillai_matrix) <- target_pillai_rownames

ae_vowels <- subset(all_vowels, VowelType %in% c("a", "e")) %>% mutate(VowelType = factor(VowelType))
ei_vowels <- subset(all_vowels, VowelType %in% c("e", "i")) %>% mutate(VowelType = factor(VowelType))
iu_vowels <- subset(all_vowels, VowelType %in% c("i", "u")) %>% mutate(VowelType = factor(VowelType))
uo_vowels <- subset(all_vowels, VowelType %in% c("u", "o")) %>% mutate(VowelType = factor(VowelType))
oa_vowels <- subset(all_vowels, VowelType %in% c("o", "a")) %>% mutate(VowelType = factor(VowelType))

edges <- list(ae_vowels, ei_vowels, iu_vowels, uo_vowels, oa_vowels)
complete_cmb_list <- list()

for (speaker in unique(aromanian$Speaker)) {
  cmb <- list()
  
  cmb[["a-e"]] <- manova(cbind(F1, F2, F3) ~ VowelType, data = subset(edges[[1]], Speaker == speaker))
  cmb[["e-i"]] <- manova(cbind(F1, F2, F3) ~ VowelType, data = subset(edges[[2]], Speaker == speaker))
  cmb[["i-u"]] <- manova(cbind(F1, F2, F3) ~ VowelType, data = subset(edges[[3]], Speaker == speaker))
  cmb[["u-o"]] <- manova(cbind(F1, F2, F3) ~ VowelType, data = subset(edges[[4]], Speaker == speaker))
  cmb[["o-a"]] <- manova(cbind(F1, F2, F3) ~ VowelType, data = subset(edges[[5]], Speaker == speaker))

  complete_cmb_list[[speaker]] <- cmb
}

other_vowel_pairs <- c("a-e", "e-i", "i-u", "u-o", "o-a")

other_pillai_results <- c()
for (speaker in unique(aromanian$Speaker)) {
  for (fct in other_vowel_pairs) {
    summary_mnv <- summary(complete_cmb_list[[speaker]][[fct]])
    other_pillai <- summary_mnv$stats[, "Pillai"]["VowelType"]
    other_pillai_results <- c(other_pillai_results, other_pillai)
  }
}

other_pillai_matrix <- matrix(other_pillai_results, nrow = length(other_vowel_pairs), ncol = length(unique(aromanian$Speaker)))

colnames(other_pillai_matrix) <- unique(aromanian$Speaker)

other_pillai_rownames <- c()
for (fct in other_vowel_pairs) {
  other_pillai_rownames <- c(other_pillai_rownames, paste(fct, "Pillai"))
}
rownames(other_pillai_matrix) <- other_pillai_rownames

means_vector <- c()
for (speaker in unique(aromanian$Speaker)) {
  five_values <- other_pillai_matrix[, speaker][1:5]
  means_vector <- c(means_vector, mean(five_values))
}
mean_pillai_matrix <- t(data.frame(means_vector))
colnames(mean_pillai_matrix) <- unique(aromanian$Speaker)


### Considering factor of vowel stress

stress_convert <- function(stress_tag) {
  if (stress_tag == "yes") {
    return("tonic")
  } else if (stress_tag == "no") {
    return("atonic")
  } else if (stress_tag == "pre") {
    return("pretonic")
  } else if (stress_tag == "post") {
    return("posttonic")
  }
}

aromanian_tonic <- subset(aromanian, Stress == "yes")
aromanian_tonic <- subset(aromanian_tonic, !(Speaker %in% c("WK_1", "WA_3", "WK_2", "WK_3", "MA_2"))) # remove raus kwim
aromanian_tonic$Speaker <- factor(aromanian_tonic$Speaker)


# Pillai scores now only on stressed occurrences

factors <- c("VowelType")
complete_cmb_list <- list()
for (speaker in unique(aromanian_tonic$Speaker)) {
  current_subset <- subset(aromanian_tonic, Speaker == speaker)
  cmb <- list()
  
  cmb[[factors[1]]] <- manova(cbind(F1, F2, F3) ~ VowelType, data = current_subset)
  
  complete_cmb_list[[speaker]] <- cmb
}

tonic_pillai_results <- c()
for (speaker in unique(aromanian_tonic$Speaker)) {
  print(speaker)
  for (fct in factors) {
    summary_mnv <- summary(complete_cmb_list[[speaker]][[fct]])
    tonic_pillai <- summary_mnv$stats[, "Pillai"][str_replace_all(fct, "\\+D\\.", "")]
    tonic_pillai_results <- c(tonic_pillai_results, tonic_pillai)
  }
}

tonic_pillai_matrix <- matrix(tonic_pillai_results, nrow = length(factors), ncol = length(unique(aromanian_tonic$Speaker)))

colnames(tonic_pillai_matrix) <- unique(aromanian_tonic$Speaker)

tonic_pillai_rownames <- c()
for (fct in factors) {
  tonic_pillai_rownames <- c(tonic_pillai_rownames, paste(fct, "Pillai"))
}
rownames(tonic_pillai_matrix) <- tonic_pillai_rownames

tonic_pillai_df <- data.frame(t(tonic_pillai_matrix))
colnames(tonic_pillai_df) <- factors
tonic_pillai_df$Speaker <- factor(rownames(tonic_pillai_df), levels = rev(rownames(tonic_pillai_df)))
rownames(tonic_pillai_df) <- NULL

tonic_speaker_exps <- expression(WT[1], WT[2], MT[1], MT[2],
                                 WA[1], WA[2], MA[1], MG[1], MS[1])


### Plotting vowel space, only 9 of 14 speakers with enough data for including stress/no-stress distinction

target_tonics <- c("/\u0259/-tonic", "/\u0268/-tonic")
target_vowels <- c("/\u0259/", "/\u0268/")
other_vowels <- c("/a/", "/e/", "/i/", "/u/", "/o/")

south_w_f1_means <- matrix(vector(), 0, 9)
colnames(south_w_f1_means) <- c(target_tonics, target_vowels, other_vowels)
south_w_f2_means <- matrix(vector(), 0, 9)
colnames(south_w_f2_means) <- c(target_tonics, target_vowels, other_vowels)

south_m_f1_means <- matrix(vector(), 0, 9)
colnames(south_m_f1_means) <- c(target_tonics, target_vowels, other_vowels)
south_m_f2_means <- matrix(vector(), 0, 9)
colnames(south_m_f2_means) <- c(target_tonics, target_vowels, other_vowels)

north_w_f1_means <- matrix(vector(), 0, 9)
colnames(north_w_f1_means) <- c(target_tonics, target_vowels, other_vowels)
north_w_f2_means <- matrix(vector(), 0, 9)
colnames(north_w_f2_means) <- c(target_tonics, target_vowels, other_vowels)

north_m_f1_means <- matrix(vector(), 0, 9)
colnames(north_m_f1_means) <- c(target_tonics, target_vowels, other_vowels)
north_m_f2_means <- matrix(vector(), 0, 9)
colnames(north_m_f2_means) <- c(target_tonics, target_vowels, other_vowels)

for (speaker in unique(aromanian_tonic$Speaker)) {
  tonic_subset <- subset(aromanian_tonic, Speaker == speaker)
  current_subset <- subset(aromanian, Speaker == speaker)
  other_subset <- subset(all_vowels, Speaker == speaker)
  f1_mean <- mean(current_subset$F1)
  f2_mean <- mean(current_subset$F2)
  dialect_and_sex <- c(as.vector(unique(current_subset$Dialect)),
                      as.vector(unique(current_subset$Sex)))
  if (all(dialect_and_sex == c("Pindian", "Female"))) {
    new_row <- c()
    new_row <- c(new_row,
                 mean(subset(tonic_subset, VowelType == 2)$F1),
                 mean(subset(tonic_subset, VowelType == 1)$F1))
    new_row <- c(new_row,
                 mean(subset(current_subset, VowelType == 2)$F1),
                 mean(subset(current_subset, VowelType == 1)$F1))
    for (vowel in other_vowels) {
      new_row <- c(new_row, mean(subset(other_subset, VowelType == substr(vowel, 2, 2))$F1))
    }
    south_w_f1_means <- rbind(south_w_f1_means, new_row)
    ###
    new_row <- c()
    new_row <- c(new_row,
                 mean(subset(tonic_subset, VowelType == 2)$F2),
                 mean(subset(tonic_subset, VowelType == 1)$F2))
    new_row <- c(new_row,
                 mean(subset(current_subset, VowelType == 2)$F2),
                 mean(subset(current_subset, VowelType == 1)$F2))
    for (vowel in other_vowels) {
      new_row <- c(new_row, mean(subset(other_subset, VowelType == substr(vowel, 2, 2))$F2))
    }
    south_w_f2_means <- rbind(south_w_f2_means, new_row)
  } else if (all(dialect_and_sex == c("Pindian", "Male"))) {
    new_row <- c()
    new_row <- c(new_row,
                 mean(subset(tonic_subset, VowelType == 2)$F1),
                 mean(subset(tonic_subset, VowelType == 1)$F1))
    new_row <- c(new_row,
                 mean(subset(current_subset, VowelType == 2)$F1),
                 mean(subset(current_subset, VowelType == 1)$F1))
    for (vowel in other_vowels) {
      new_row <- c(new_row, mean(subset(other_subset, VowelType == substr(vowel, 2, 2))$F1))
    }
    south_m_f1_means <- rbind(south_m_f1_means, new_row)
    ###
    new_row <- c()
    new_row <- c(new_row,
                 mean(subset(tonic_subset, VowelType == 2)$F2),
                 mean(subset(tonic_subset, VowelType == 1)$F2))
    new_row <- c(new_row,
                 mean(subset(current_subset, VowelType == 2)$F2),
                 mean(subset(current_subset, VowelType == 1)$F2))
    for (vowel in other_vowels) {
      new_row <- c(new_row, mean(subset(other_subset, VowelType == substr(vowel, 2, 2))$F2))
    }
    south_m_f2_means <- rbind(south_m_f2_means, new_row)
  } else if (all(dialect_and_sex == c("Farsherot", "Female"))) {
    new_row <- c()
    new_row <- c(new_row,
                 mean(subset(tonic_subset, VowelType == 2)$F1),
                 mean(subset(tonic_subset, VowelType == 1)$F1))
    new_row <- c(new_row,
                 mean(subset(current_subset, VowelType == 2)$F1),
                 mean(subset(current_subset, VowelType == 1)$F1))
    for (vowel in other_vowels) {
      new_row <- c(new_row, mean(subset(other_subset, VowelType == substr(vowel, 2, 2))$F1))
    }
    north_w_f1_means <- rbind(north_w_f1_means, new_row)
    ###
    new_row <- c()
    new_row <- c(new_row,
                 mean(subset(tonic_subset, VowelType == 2)$F2),
                 mean(subset(tonic_subset, VowelType == 1)$F2))
    new_row <- c(new_row,
                 mean(subset(current_subset, VowelType == 2)$F2),
                 mean(subset(current_subset, VowelType == 1)$F2))
    for (vowel in other_vowels) {
      new_row <- c(new_row, mean(subset(other_subset, VowelType == substr(vowel, 2, 2))$F2))
    }
    north_w_f2_means <- rbind(north_w_f2_means, new_row)
  } else if (all(dialect_and_sex == c("Farsherot", "Male"))) {
    new_row <- c()
    new_row <- c(new_row,
                 mean(subset(tonic_subset, VowelType == 2)$F1),
                 mean(subset(tonic_subset, VowelType == 1)$F1))
    new_row <- c(new_row,
                 mean(subset(current_subset, VowelType == 2)$F1),
                 mean(subset(current_subset, VowelType == 1)$F1))
    for (vowel in other_vowels) {
      new_row <- c(new_row, mean(subset(other_subset, VowelType == substr(vowel, 2, 2))$F1))
    }
    north_m_f1_means <- rbind(north_m_f1_means, new_row)
    ###
    new_row <- c()
    new_row <- c(new_row,
                 mean(subset(tonic_subset, VowelType == 2)$F2),
                 mean(subset(tonic_subset, VowelType == 1)$F2))
    new_row <- c(new_row,
                 mean(subset(current_subset, VowelType == 2)$F2),
                 mean(subset(current_subset, VowelType == 1)$F2))
    for (vowel in other_vowels) {
      new_row <- c(new_row, mean(subset(other_subset, VowelType == substr(vowel, 2, 2))$F2))
    }
    north_m_f2_means <- rbind(north_m_f2_means, new_row)
  }
}

south_info <- matrix(vector(), 0, 4)
colnames(south_info) <- c("Vowel", "Sex", "F1", "F2")

north_info <- matrix(vector(), 0, 4)
colnames(north_info) <- c("Vowel", "Sex", "F1", "F2")

for (vowel in c(target_tonics, target_vowels, other_vowels)) {
  south_info <- rbind(south_info, c(vowel, "Female", mean(south_w_f1_means[, vowel]), mean(south_w_f2_means[, vowel])))
  south_info <- rbind(south_info, c(vowel, "Male", mean(south_m_f1_means[, vowel]), mean(south_m_f2_means[, vowel])))
}

for (vowel in c(target_tonics, target_vowels, other_vowels)) {
  north_info <- rbind(north_info, c(vowel, "Female", mean(north_w_f1_means[, vowel]), mean(north_w_f2_means[, vowel])))
  north_info <- rbind(north_info, c(vowel, "Male", mean(north_m_f1_means[, vowel]), mean(north_m_f2_means[, vowel])))
}

south_info <- data.frame(south_info)
north_info <- data.frame(north_info)
south_info <- mutate(south_info, Vowel = as.factor(Vowel), Sex = as.factor(Sex), F1 = as.numeric(F1), F2 = as.numeric(F2))
north_info <- mutate(north_info, Vowel = as.factor(Vowel), Sex = as.factor(Sex), F1 = as.numeric(F1), F2 = as.numeric(F2))

customPalette <- c("black", "#37B98D")

south_w_plot <- ggplot(subset(south_info, Sex == "Female"), aes(x = F2, y = F1, color = ifelse(!is.na(str_match(Vowel, "tonic")), "yes", "no"))) +
  geom_text(aes(label = substring(Vowel, 1, 3)), size = 3.5) +
  scale_x_reverse() + scale_y_reverse() +
  labs(x = "F2 (Bark)", y = "F1 (Bark)") +
  ggtitle(expression(paste(WT[1], ", ", WT[2], sep="")),
          "Avg. Pindian female vowel space") +
  scale_color_manual(values=customPalette) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  theme(legend.position = "none")

south_m_plot <- ggplot(subset(south_info, Sex == "Male"), aes(x = F2, y = F1, color = ifelse(!is.na(str_match(Vowel, "tonic")), "yes", "no"))) +
  geom_text(aes(label = substring(Vowel, 1, 3)), size = 3.5) +
  scale_x_reverse() + scale_y_reverse() +
  labs(x = "F2 (Bark)", y = "F1 (Bark)") +
  ggtitle(expression(paste(MT[1], ", ", MT[2], sep="")),
          "Avg. Pindian male vowel space") +
  scale_color_manual(values=customPalette) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  theme(legend.position = "none")

north_w_plot <- ggplot(subset(north_info, Sex == "Female"), aes(x = F2, y = F1, color = ifelse(!is.na(str_match(Vowel, "tonic")), "yes", "no"))) +
  geom_text(aes(label = substring(Vowel, 1, 3)), size = 3.5) +
  scale_x_reverse() + scale_y_reverse() +
  labs(x = "F2 (Bark)", y = "F1 (Bark)", color = "Stressed only?") +
  ggtitle(expression(paste(WA[1], ", ", WA[2], sep="")),
          "Avg. Farsherot female vowel space") +
  scale_color_manual(values=customPalette) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  theme(legend.position = c(.74, .36),
        legend.justification = c("center", "center"),
        legend.box.just = "right",
        legend.margin = margin(12, 18, 12, 18),
        legend.background = element_rect(fill=alpha('gray90', 1))) +
  guides(color = guide_legend(reverse = TRUE)) 

north_m_plot <- ggplot(subset(north_info, Sex == "Male"), aes(x = F2, y = F1, color = ifelse(!is.na(str_match(Vowel, "tonic")), "yes", "no"))) +
  geom_text(aes(label = substring(Vowel, 1, 3)), size = 3.5) +
  scale_x_reverse() + scale_y_reverse() +
  labs(x = "F2 (Bark)", y = "F1 (Bark)") +
  ggtitle(expression(paste(MA[1], ", ", MG[1], ", ", MS[1], sep="")),
          "Avg. Farsherot male vowel space") +
  scale_color_manual(values=customPalette) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  theme(legend.position = "none")

grid.arrange(rbind(cbind(ggplotGrob(south_w_plot), ggplotGrob(north_w_plot), size="last"),
                   cbind(ggplotGrob(south_m_plot), ggplotGrob(north_m_plot), size="last")))
