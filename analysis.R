library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(emuR)
library(visdat)
library(stringr)
library(car)
library(rstatix)
library(gridExtra)
library(Cairo)
library(heplots)

# Some functions:

merge_contexts_common <- function(phone) {
  if (phone %in% c("v", "f")) {
    return("v/f")
  } else if (phone %in% c("D", "T")) {
    return("D/T")
  } else if (phone %in% c("G", "R", "x")) {
    return("G/R/x")
  } else if (phone %in% c("b", "p")) {
    return("b/p")
  } else {
    return(phone)
  }
}

merge_contexts_before <- function(phone) {
  if (phone %in% c("z", "s", "dz", "ts")) {
    return("z/s/dz/ts_")
  } else if (phone %in% c("g", "k")) {
    return("g/k_") 
  } else if (phone %in% c("Z", "S", "dZ", "tS")) {
    return("Z/S/dZ/tS_")
  } else if (phone %in% c("d", "t")) {
    return("d/t_")
  } else if (phone %in% c("j", "J|", "c", "C")) {
    return("j/J|/c/C_")
  } else {
    return(paste(phone, "_", sep=""))
  }
}

merge_contexts_after <- function(phone) {
  if (phone %in% c("d", "t", "dz", "ts", "dZ", "tS")) {
    return("_d/t/dz/ts/dZ/tS")
  } else if (phone %in% c("J|", "g", "c", "k")) {
    return("_J|/g/c/k")
  } else if (phone %in% c("z", "s")) {
    return("_z/s")
  } else if (phone %in% c("j", "C")) {
    return("_j/C")
  } else if (phone %in% c("Z", "S")) {
    return("_Z/S")
  } else {
    return(paste("_", phone, sep=""))
  }
}

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

cbPalette <- c("#56B4E9", "#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00", "chartreuse3", "#999999")

# Initial preparation of datasets:

aromanian <- read.csv('C:/Users/Michael/OneDrive/Stuff/wise2021/ba/aromanian_targets.csv')
all_vowels <- read.csv('C:/Users/Michael/OneDrive/Stuff/wise2021/ba/aromanian_all_raw.csv')

aromanian <- separate(aromanian, col = CoArtContext, into = c("LeftCoArt", "RightCoArt"), sep = "_")
all_vowels <- separate(all_vowels, col = CoArtContext, into = c("LeftCoArt", "RightCoArt"), sep = "_")

print(union(unique(aromanian$LeftCoArt), unique(aromanian$RightCoArt)))
print(length(union(unique(aromanian$LeftCoArt), unique(aromanian$RightCoArt))))
print(union(unique(all_vowels$LeftCoArt), unique(all_vowels$RightCoArt)))
print(length(union(unique(all_vowels$LeftCoArt), unique(all_vowels$RightCoArt))))

aromanian$LeftCoArt <- sapply(sapply(aromanian$LeftCoArt, merge_contexts_common), merge_contexts_before)
all_vowels$LeftCoArt <- sapply(sapply(all_vowels$LeftCoArt, merge_contexts_common), merge_contexts_before)

aromanian$RightCoArt <- sapply(sapply(aromanian$RightCoArt, merge_contexts_common), merge_contexts_after)
all_vowels$RightCoArt <- sapply(sapply(all_vowels$RightCoArt, merge_contexts_common), merge_contexts_after)

aromanian$Sex <- ifelse(substring(aromanian$Speaker, 1, 1) == "W", "Female", "Male")
all_vowels$Sex <- ifelse(substring(all_vowels$Speaker, 1, 1) == "W", "Female", "Male")

aromanian$InitialID <- sapply(sapply(seq(nrow(aromanian)), toString), function(x) paste0("T", str_pad(x, 4, pad = "0")))
all_vowels$InitialID <- sapply(sapply(seq(nrow(all_vowels)), toString), function(x) paste0("V", str_pad(x, 5, pad = "0")))

aromanian <- mutate(aromanian, VowelType = as.factor(VowelType),
                    Region = as.factor(Region), Speaker = as.factor(Speaker),
                    POS = as.factor(POS), Stress = as.factor(Stress),
                    LeftCoArt = as.factor(LeftCoArt),
                    RightCoArt = as.factor(RightCoArt))
all_vowels <- mutate(all_vowels, VowelType = as.factor(substring(str_extract(VowelInWord, "\\(.+\\)"), 2, 2)),
                     Region = as.factor(Region), Speaker = as.factor(Speaker),
                     LeftCoArt = as.factor(LeftCoArt),
                     RightCoArt = as.factor(RightCoArt))

glimpse(aromanian) # there you show in paper
glimpse(all_vowels) # there you show in paper

aromanian <- mutate(aromanian, F1 = bark(F1), F2 = bark(F2), F3 = bark(F3))
all_vowels <- mutate(all_vowels, F1 = bark(F1), F2 = bark(F2), F3 = bark(F3))

### here create some plots TODO:
plots <- list()
i = 0
for (current_speaker in c("MT_1", "WA_1")) {
  i <- i + 1
  current_subset <- subset(aromanian, Speaker == current_speaker & VowelType == "2")
  formants_info <- rbind(cbind("F1", current_subset$F1),
                         cbind("F2", current_subset$F2),
                         cbind("F3", current_subset$F3))
  formants_info <- data.frame(formants_info)
  colnames(formants_info) <- c("Formant", "Bark")
  formants_info <- mutate(formants_info, Formant = as.factor(Formant),
                          Bark = as.numeric(Bark))
  current_plot <- ggplot(formants_info, aes(x = Formant, y = Bark)) +
                    geom_boxplot() +
                    ylim(2, 18) +
                    theme_classic() +
                    theme(plot.title = element_text(hjust = 0.5),
                          plot.subtitle = element_text(hjust = 0.5))
  if (current_speaker == "MT_1") {
    current_plot <- current_plot + ggtitle(expression(paste("Speaker ", MT[1], sep="")),
                                           "Realizations of /\u0259/")
  } else {
    current_plot <- current_plot + ggtitle(expression(paste("Speaker ", WA[1], sep="")),
                                           "Realizations of /\u0259/")
  }
  plots[[i]] <- current_plot
}

grid.arrange(cbind(ggplotGrob(plots[[1]]), ggplotGrob(plots[[2]]), size="last"))

###
already_filtered <- FALSE
###
# MD outlier exclusion:
plots <- list() ### for TODO below
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
      
      ### create some plots TODO
      if (speaker %in% c("MT_1", "WA_1") & vowel == "2") {
        if (speaker == "MT_1") {
          View(subset(current_subset, Speaker == "MT_1"))
        }
        other_i <- other_i + 1
        current_plot <- ggplot(current_subset, aes(MDBySpeaker)) +
                          geom_histogram(bins = ifelse(speaker == "MT_1", 35, 15), na.rm = T) +
                          geom_vline(aes(xintercept=7.82),
                                     color="red", linetype="dashed", size=0.5) +
                          xlim(0, ifelse(speaker == "MT_1", 35, 15)) +
                          ylim(0, ifelse(speaker == "MT_1", 55, 35)) +
                          labs(x = "Squared Mahalanobis distance", y = "Count") +
                          theme_classic() +
                          theme(plot.title = element_text(hjust = 0.5),
                                plot.subtitle = element_text(hjust = 0.5))
        if (speaker == "MT_1") {
          current_plot <- current_plot + ggtitle(expression(paste("Speaker ", MT[1], sep="")),
                                                 "Realizations of /\u0259/")
        } else {
          current_plot <- current_plot + ggtitle(expression(paste("Speaker ", WA[1], sep="")),
                                                 "Realizations of /\u0259/")
        }
        plots[[other_i]] <- current_plot
      }
      ###
      
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
  
  # Show the ones that are filtered out for targets:
  print(filter(aromanian, IsExtremeOut == TRUE)$InitialID)
  
  # Show the ones that are filtered out for all:
  print(filter(all_vowels, IsExtremeOut == TRUE)$InitialID)
  
  ### create some plots TODO
  plots_new <- list()
  i = 0
  for (current_speaker in c("MT_1", "WA_1")) {
    i <- i + 1
    current_subset <- subset(aromanian, Speaker == current_speaker & VowelType == "2")
    current_plot <- ggplot(current_subset, aes(x = F2, y = F1, color = IsExtremeOut)) +
                      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
                      geom_point(alpha=0.6) +
                      xlim(16, 8) + ylim(10, 2) +
                      labs(x = "F2 (Bark)", y = "F1 (Bark)") +
                      theme_classic() +
                      theme(legend.position = "none") +
                      theme(plot.title = element_text(hjust = 0.5),
                            plot.subtitle = element_text(hjust = 0.5))
    if (current_speaker == "MT_1") {
      current_plot <- current_plot + ggtitle(expression(paste("Speaker ", MT[1], sep="")),
                                             "Realizations of /\u0259/")
    } else {
      current_plot <- current_plot + ggtitle(expression(paste("Speaker ", WA[1], sep="")),
                                             "Realizations of /\u0259/")
    }
    plots_new[[i * 2 - 1]] <- current_plot
    current_plot <- ggplot(current_subset, aes(x = F2, y = F3, color = ifelse(IsExtremeOut, "yes", "no"))) +
      scale_color_manual(values = c("yes" = "red", "no" = "black")) +
      geom_point(alpha=0.6) +
      xlim(16, 8) + ylim(20, 12) +
      labs(x = "F2 (Bark)", y = "F3 (Bark)") +
      theme_classic() +
      guides(color = guide_legend(reverse = TRUE))
    if (i == 2) {
      current_plot <- current_plot +
                        theme(legend.position = c(.9, .4),
                              legend.justification = c("right", "top"),
                              legend.box.just = "right",
                              legend.margin = margin(12, 18, 12, 18),
                              legend.background = element_rect(fill=alpha('gray90', 1))) +
                        labs(color = 'Outlier?')
    } else {
      current_plot <- current_plot + theme(legend.position = "none")
    }
    plots_new[[i * 2]] <- current_plot
    print(current_plot)
  }
  ###

  aromanian <- subset(aromanian, IsExtremeOut == FALSE)
  aromanian <- subset(aromanian, select = -c(IsExtremeOut, IsModerateOut))
  
  all_vowels <- subset(all_vowels, IsExtremeOut == FALSE)
  all_vowels <- subset(all_vowels, select = -c(IsExtremeOut, IsModerateOut))
  
  already_filtered = TRUE
}
grid.arrange(cbind(ggplotGrob(plots[[1]]), ggplotGrob(plots[[2]]), size="last"))

grid.arrange(rbind(cbind(ggplotGrob(plots_new[[1]]), ggplotGrob(plots_new[[3]]), size="last"),
                   cbind(ggplotGrob(plots_new[[2]]), ggplotGrob(plots_new[[4]]), size="last")))


# Excluding pre- and post-vocalic contexts for all:

glimpse(aromanian)

glimpse(all_vowels)

### manova stuff
factors <- c("POS", "LeftCoArt", "RightCoArt", "Stress", "VowelType",
             "POS+D.", "LeftCoArt+D.", "RightCoArt+D.", "Stress+D.", "VowelType+D.")

complete_cmb_list <- list()
for (speaker in unique(aromanian$Speaker)) {
  print(speaker)
  current_subset <- subset(aromanian, Speaker == speaker)
  cmb <- list()

  ## like this:
  cmb[[factors[1]]] <- manova(cbind(F1, F2, F3) ~ POS, data = current_subset)
  cmb[[factors[2]]] <- manova(cbind(F1, F2, F3) ~ LeftCoArt, data = current_subset)
  cmb[[factors[3]]] <- manova(cbind(F1, F2, F3) ~ RightCoArt, data = current_subset)
  cmb[[factors[4]]] <- manova(cbind(F1, F2, F3) ~ Stress, data = current_subset)
  cmb[[factors[5]]] <- manova(cbind(F1, F2, F3) ~ VowelType, data = current_subset)
  ## controling for duration:
  cmb[[factors[6]]] <- manova(cbind(F1, F2, F3) ~ POS + Duration, data = current_subset)
  cmb[[factors[7]]] <- manova(cbind(F1, F2, F3) ~ LeftCoArt + Duration, data = current_subset)
  cmb[[factors[8]]] <- manova(cbind(F1, F2, F3) ~ RightCoArt + Duration, data = current_subset)
  cmb[[factors[9]]] <- manova(cbind(F1, F2, F3) ~ Stress + Duration, data = current_subset)
  cmb[[factors[10]]] <- manova(cbind(F1, F2, F3) ~ VowelType + Duration, data = current_subset)

  complete_cmb_list[[speaker]] <- cmb
}

pillai_results <- c()
p_results <- c()
for (speaker in unique(aromanian$Speaker)) {
  for (fct in factors) {
    summary_mnv <- summary(complete_cmb_list[[speaker]][[fct]])
    pillai <- summary_mnv$stats[, "Pillai"][str_replace_all(fct, "\\+D\\.", "")]
    p <- summary_mnv$stats[, "Pr(>F)"][str_replace_all(fct, "\\+D\\.", "")]
    pillai_results <- c(pillai_results, pillai)
    p_results <- c(p_results, p)
  }
}

pillai_matrix <- matrix(pillai_results, nrow = length(factors), ncol = length(unique(aromanian$Speaker)))
p_matrix <- matrix(p_results, nrow = length(factors), ncol = length(unique(aromanian$Speaker)))

colnames(pillai_matrix) <- unique(aromanian$Speaker)
colnames(p_matrix) <- unique(aromanian$Speaker)

pillai_rownames <- c()
p_rownames <- c() 
for (fct in factors) {
  pillai_rownames <- c(pillai_rownames, paste(fct, "Pillai"))
  p_rownames <- c(p_rownames, paste(fct, "p-value"))
}
rownames(pillai_matrix) <- pillai_rownames
rownames(p_matrix) <- p_rownames

significance_matrix <- p_matrix < 0.05

pillai_matrix_signif <- matrix(nrow = nrow(p_matrix), ncol = ncol(p_matrix))
p_matrix_signif <- matrix(nrow = nrow(p_matrix), ncol = ncol(p_matrix))

colnames(pillai_matrix_signif) <- unique(aromanian$Speaker)
colnames(p_matrix_signif) <- unique(aromanian$Speaker)
rownames(pillai_matrix_signif) <- pillai_rownames
rownames(p_matrix_signif) <- p_rownames

for (i in 1:nrow(p_matrix)) {
  for (j in 1:ncol(p_matrix)) {
    if (!is.na(significance_matrix[i, j]) & significance_matrix[i, j]) {
      pillai_matrix_signif[i, j] = pillai_matrix[i, j]
      p_matrix_signif[i, j] = p_matrix[i, j]
    }
  }
}

write.csv(data.frame(pillai_matrix), "pillai_matrix.csv")

write.csv(data.frame(p_matrix), "p_matrix.csv")

write.csv(data.frame(pillai_matrix_signif), "pillai_matrix_signif.csv")
###

'
current_subset <- subset(aromanian, Speaker == "WK_1")
ggplot(current_subset, aes(x=F2, y=F1, color=VowelType)) +
        geom_text(aes(label=VowelInWord)) +
        scale_x_reverse() + scale_y_reverse()
'

pillai_df <- data.frame(t(pillai_matrix))
colnames(pillai_df) <- factors
pillai_df$Speaker <- factor(rownames(pillai_df), levels = rev(rownames(pillai_df)))
rownames(pillai_df) <- NULL
pillai_likethis_df <- subset(pillai_df, select = c(factors[1:5], "Speaker"))
pillai_durcontrol_df <- subset(pillai_df, select = c(factors[6:10], "Speaker"))
colnames(pillai_durcontrol_df) <- c(factors[1:5], "Speaker")
pillai_likethis_df$DurControl <- FALSE
pillai_durcontrol_df$DurControl <- TRUE
pillai_df <- rbind(pillai_likethis_df, pillai_durcontrol_df)

speaker_exps <- expression(WT[1], WT[2], WK[1], WK[2], WK[3], MT[1], MT[2],
                           WA[1], WA[2], WA[3], MA[1], MA[2], MG[1], MS[1])

#cbPalette <- c("#56B4E9", "#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00", "chartreuse3", "#999999")
customPalette <- c("#E69F00", "#999999")


ggplot(pillai_df, aes(x = Speaker, y = VowelType, fill = ifelse(DurControl, "yes", "no"))) +
  geom_bar(position="dodge", stat="identity", width = 0.8) +
  labs(x = "Speaker", y = "Pillai score between /\u0259/ and /\u0268/") +
  ylim(0, 1) +
  coord_flip() +
  scale_x_discrete(breaks = unique(aromanian$Speaker),
                   labels = speaker_exps) +
  scale_fill_manual(values=customPalette) +
  theme_classic() +
  theme(legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(12, 18, 12, 18),
        legend.background = element_rect(fill=alpha('gray90', 1))) +
  labs(fill = 'Controlled for duration?') +
  guides(fill = guide_legend(reverse = TRUE)) 

print("Mean pillai South:")
print(mean(pillai_likethis_df$VowelType[1:7]))

print("Mean pillai North:")
print(mean(pillai_likethis_df$VowelType[8:14]))

print("Mean pillai South without WK_2:")
print(mean(c(pillai_likethis_df$VowelType[1:3], pillai_likethis_df$VowelType[5:7])))

print("Mean pillai North without WA_3:")
print(mean(c(pillai_likethis_df$VowelType[8:9], pillai_likethis_df$VowelType[11:14])))

rem_lengths <- c(150, 82, 32, 16, 36, 176, 50,
                 152, 116, 15, 67, 41, 369, 145)
cor(rem_lengths, pillai_likethis_df$VowelType)

# pentagram solution:
ae_vowels <- subset(all_vowels, VowelType %in% c("a", "e")) %>% mutate(VowelType = factor(VowelType))
ei_vowels <- subset(all_vowels, VowelType %in% c("e", "i")) %>% mutate(VowelType = factor(VowelType))
iu_vowels <- subset(all_vowels, VowelType %in% c("i", "u")) %>% mutate(VowelType = factor(VowelType))
uo_vowels <- subset(all_vowels, VowelType %in% c("u", "o")) %>% mutate(VowelType = factor(VowelType))
oa_vowels <- subset(all_vowels, VowelType %in% c("o", "a")) %>% mutate(VowelType = factor(VowelType))

### pentagram manova stuff ###
edges <- list(ae_vowels, ei_vowels, iu_vowels, uo_vowels, oa_vowels)
complete_cmb_list <- list()

for (speaker in unique(aromanian$Speaker)) {
  #print(speaker)
  cmb <- list()
  # two-way
  ## like this:
  cmb[["a-e"]] <- manova(cbind(F1, F2, F3) ~ VowelType, data = subset(edges[[1]], Speaker == speaker))
  cmb[["e-i"]] <- manova(cbind(F1, F2, F3) ~ VowelType, data = subset(edges[[2]], Speaker == speaker))
  cmb[["i-u"]] <- manova(cbind(F1, F2, F3) ~ VowelType, data = subset(edges[[3]], Speaker == speaker))
  cmb[["u-o"]] <- manova(cbind(F1, F2, F3) ~ VowelType, data = subset(edges[[4]], Speaker == speaker))
  cmb[["o-a"]] <- manova(cbind(F1, F2, F3) ~ VowelType, data = subset(edges[[5]], Speaker == speaker))
  ## controling for duration:
  cmb[["a-e+D."]] <- manova(cbind(F1, F2, F3) ~ VowelType + Duration, data = subset(edges[[1]], Speaker == speaker))
  cmb[["e-i+D."]] <- manova(cbind(F1, F2, F3) ~ VowelType + Duration, data = subset(edges[[2]], Speaker == speaker))
  cmb[["i-u+D."]] <- manova(cbind(F1, F2, F3) ~ VowelType + Duration, data = subset(edges[[3]], Speaker == speaker))
  cmb[["u-o+D."]] <- manova(cbind(F1, F2, F3) ~ VowelType + Duration, data = subset(edges[[4]], Speaker == speaker))
  cmb[["o-a+D."]] <- manova(cbind(F1, F2, F3) ~ VowelType + Duration, data = subset(edges[[5]], Speaker == speaker))
  
  complete_cmb_list[[speaker]] <- cmb
}

things <- c("a-e", "e-i", "i-u", "u-o", "o-a",
            "a-e+D.", "e-i+D.", "i-u+D.", "u-o+D.", "o-a+D.")

pent_pillai_results <- c()
pent_p_results <- c()
for (speaker in unique(aromanian$Speaker)) {
  for (fct in things) {
    summary_mnv <- summary(complete_cmb_list[[speaker]][[fct]])
    pent_pillai <- summary_mnv$stats[, "Pillai"]["VowelType"]
    pent_p <- summary_mnv$stats[, "Pr(>F)"]["VowelType"]
    pent_pillai_results <- c(pent_pillai_results, pent_pillai)
    pent_p_results <- c(pent_p_results, pent_p)
  }
}

pent_pillai_matrix <- matrix(pent_pillai_results, nrow = length(things), ncol = length(unique(aromanian$Speaker)))
pent_p_matrix <- matrix(pent_p_results, nrow = length(things), ncol = length(unique(aromanian$Speaker)))

colnames(pent_pillai_matrix) <- unique(aromanian$Speaker)
colnames(pent_p_matrix) <- unique(aromanian$Speaker)

pent_pillai_rownames <- c()
pent_p_rownames <- c() 
for (fct in things) {
  pent_pillai_rownames <- c(pent_pillai_rownames, paste(fct, "Pillai"))
  pent_p_rownames <- c(pent_p_rownames, paste(fct, "p-value"))
}
rownames(pent_pillai_matrix) <- pent_pillai_rownames
rownames(pent_p_matrix) <- pent_p_rownames

pent_significance_matrix <- pent_p_matrix < 0.05

pent_pillai_matrix_signif <- matrix(nrow = nrow(pent_p_matrix), ncol = ncol(pent_p_matrix))
pent_p_matrix_signif <- matrix(nrow = nrow(pent_p_matrix), ncol = ncol(pent_p_matrix))

colnames(pent_pillai_matrix_signif) <- unique(aromanian$Speaker)
colnames(pent_p_matrix_signif) <- unique(aromanian$Speaker)
rownames(pent_pillai_matrix_signif) <- pillai_rownames
rownames(pent_p_matrix_signif) <- p_rownames

for (i in 1:nrow(pent_p_matrix)) {
  for (j in 1:ncol(pent_p_matrix)) {
    if (!is.na(pent_significance_matrix[i, j]) & pent_significance_matrix[i, j]) {
      pent_pillai_matrix_signif[i, j] = pent_pillai_matrix[i, j]
      pent_p_matrix_signif[i, j] = pent_p_matrix[i, j]
    }
  }
}

write.csv(data.frame(pent_pillai_matrix), "pent_pillai_matrix.csv")

write.csv(data.frame(pent_p_matrix), "pent_p_matrix.csv")

write.csv(data.frame(pent_pillai_matrix_signif), "pent_pillai_matrix_signif.csv")

means_vector <- c()
for (speaker in unique(aromanian$Speaker)) {
  five_values <- pent_pillai_matrix[, speaker][1:5]
  means_vector <- c(means_vector, mean(five_values))
}
one_row_df <- t(data.frame(means_vector))
colnames(one_row_df) <- unique(aromanian$Speaker)


###
###
###

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

#cbPalette <- c("#56B4E9", "#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00", "chartreuse3", "#999999")
customPalette <- c("#999999", "#D55E00", "chartreuse3", "#009E73")
wa3_vowel_plot <- ggplot(subset(aromanian, Speaker == "WA_3"), aes(x = F2, y = F1, color = ifelse(VowelType == 2, "/\u0259/", "/\u0268/"))) +
                    geom_text(aes(label = VowelInWord), size = 3) +
                    xlim(15, 10) + ylim(6.5, 4) +
                    labs(x = "F2 (Bark)", y = "F1 (Bark)", color = "Vowel") +
                    ggtitle(expression(paste("Speaker ", WA[3], sep="")),
                            "Realizations of /\u0259/ and /\u0268/") +
                    scale_colour_manual(values=cbPalette) +
                    theme_classic() +
                    theme(plot.title = element_text(hjust = 0.5),
                          plot.subtitle = element_text(hjust = 0.5)) +
                    theme(legend.position = c(.8, .3),
                          legend.justification = c("center", "center"),
                          legend.box.just = "right",
                          legend.margin = margin(12, 18, 12, 18),
                          legend.background = element_rect(fill=alpha('gray90', 1)))

wa3_stress_plot <- ggplot(subset(aromanian, Speaker == "WA_3"), aes(x = F2, y = F1, color = Stress)) +
                     geom_text(aes(label = sapply(Stress, stress_convert)), size = 3) +
                     xlim(15, 10) + ylim(6.5, 4) +
                     labs(x = "F2 (Bark)", y = "F1 (Bark)", color = "Stress") +
                     ggtitle(expression(paste("Speaker ", WA[3], sep="")),
                             "Realizations of /\u0259/ and /\u0268/") +
                     scale_colour_manual(values=customPalette) +
                     theme_classic() +
                     theme(plot.title = element_text(hjust = 0.5),
                           plot.subtitle = element_text(hjust = 0.5)) +
                     theme(legend.position = "none")

grid.arrange(cbind(ggplotGrob(wa3_vowel_plot), ggplotGrob(wa3_stress_plot), size="last"))


table(aromanian$VowelType, aromanian$Stress)

# tonic only:
aromanian_tonic <- subset(aromanian, Stress == "yes")
aromanian_tonic <- subset(aromanian_tonic, !(Speaker %in% c("WK_1", "WA_3", "WK_2", "WK_3", "MA_2"))) # remove raus kwim
aromanian_tonic$Speaker <- factor(aromanian_tonic$Speaker)

### manova stuff tonic only:
factors <- c("VowelType", "VowelType+D.")
complete_cmb_list <- list()
for (speaker in unique(aromanian_tonic$Speaker)) {
  print(speaker)
  current_subset <- subset(aromanian_tonic, Speaker == speaker)
  cmb <- list()
  
  cmb[[factors[1]]] <- manova(cbind(F1, F2, F3) ~ VowelType, data = current_subset)
  cmb[[factors[2]]] <- manova(cbind(F1, F2, F3) ~ VowelType + Duration, data = current_subset)
  
  complete_cmb_list[[speaker]] <- cmb
}

tonic_pillai_results <- c()
tonic_p_results <- c()
for (speaker in unique(aromanian_tonic$Speaker)) {
  print(speaker)
  for (fct in factors) {
    summary_mnv <- summary(complete_cmb_list[[speaker]][[fct]])
    tonic_pillai <- summary_mnv$stats[, "Pillai"][str_replace_all(fct, "\\+D\\.", "")]
    tonic_p <- summary_mnv$stats[, "Pr(>F)"][str_replace_all(fct, "\\+D\\.", "")]
    tonic_pillai_results <- c(tonic_pillai_results, tonic_pillai)
    tonic_p_results <- c(tonic_p_results, tonic_p)
  }
}

tonic_pillai_matrix <- matrix(tonic_pillai_results, nrow = length(factors), ncol = length(unique(aromanian_tonic$Speaker)))
tonic_p_matrix <- matrix(tonic_p_results, nrow = length(factors), ncol = length(unique(aromanian_tonic$Speaker)))

colnames(tonic_pillai_matrix) <- unique(aromanian_tonic$Speaker)
colnames(tonic_p_matrix) <- unique(aromanian_tonic$Speaker)

tonic_pillai_rownames <- c()
tonic_p_rownames <- c() 
for (fct in factors) {
  tonic_pillai_rownames <- c(tonic_pillai_rownames, paste(fct, "Pillai"))
  tonic_p_rownames <- c(tonic_p_rownames, paste(fct, "p-value"))
}
rownames(tonic_pillai_matrix) <- tonic_pillai_rownames
rownames(tonic_p_matrix) <- tonic_p_rownames

tonic_significance_matrix <- tonic_p_matrix < 0.05

tonic_pillai_matrix_signif <- matrix(nrow = nrow(tonic_p_matrix), ncol = ncol(tonic_p_matrix))
tonic_p_matrix_signif <- matrix(nrow = nrow(tonic_p_matrix), ncol = ncol(tonic_p_matrix))

colnames(tonic_pillai_matrix_signif) <- unique(aromanian_tonic$Speaker)
colnames(tonic_p_matrix_signif) <- unique(aromanian_tonic$Speaker)
rownames(tonic_pillai_matrix_signif) <- tonic_pillai_rownames
rownames(tonic_p_matrix_signif) <- tonic_p_rownames

for (i in 1:nrow(tonic_p_matrix)) {
  for (j in 1:ncol(tonic_p_matrix)) {
    if (!is.na(tonic_significance_matrix[i, j]) & tonic_significance_matrix[i, j]) {
      tonic_pillai_matrix_signif[i, j] = tonic_pillai_matrix[i, j]
      tonic_p_matrix_signif[i, j] = tonic_p_matrix[i, j]
    }
  }
}

write.csv(data.frame(tonic_pillai_matrix), "tonic_pillai_matrix.csv")

write.csv(data.frame(tonic_p_matrix), "tonic_p_matrix.csv")

write.csv(data.frame(tonic_pillai_matrix_signif), "tonic_pillai_matrix_signif.csv")

###

#copy-pasted geom_bar stuff:

tonic_pillai_df <- data.frame(t(tonic_pillai_matrix))
colnames(tonic_pillai_df) <- factors
tonic_pillai_df$Speaker <- factor(rownames(tonic_pillai_df), levels = rev(rownames(tonic_pillai_df)))
rownames(tonic_pillai_df) <- NULL
tonic_pillai_likethis_df <- subset(tonic_pillai_df, select = c(factors[1], "Speaker"))
general_pillai_likethis_df <- subset(pillai_likethis_df, !(Speaker %in% c("WK_1", "WA_3", "WK_2", "WK_3", "MA_2")),
                                     select = c("VowelType", "Speaker"))
rownames(general_pillai_likethis_df) <- NULL
tonic_pillai_likethis_df$TonicOnly <- TRUE
general_pillai_likethis_df$TonicOnly <- FALSE
tonic_pillai_df <- rbind(tonic_pillai_likethis_df, general_pillai_likethis_df)

tonic_speaker_exps <- expression(WT[1], WT[2], MT[1], MT[2],
                                 WA[1], WA[2], MA[1], MG[1], MS[1])

#cbPalette <- c("#56B4E9", "#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00", "chartreuse3", "#999999")
customPalette <- c("#E69F00", "#009E73")

ggplot(tonic_pillai_df, aes(x = Speaker, y = VowelType, fill = ifelse(TonicOnly, "yes", "no"))) +
  geom_bar(position="dodge", stat="identity", width = 0.8) +
  labs(x = "Speaker", y = "Pillai score between /\u0259/ and /\u0268/") +
  ylim(0, 1) +
  coord_flip() +
  scale_x_discrete(breaks = unique(aromanian$Speaker),
                   labels = speaker_exps) +
  scale_fill_manual(values=customPalette) +
  theme_classic() +
  theme(legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(12, 18, 12, 18),
        legend.background = element_rect(fill=alpha('gray90', 1))) +
  labs(fill = 'Tonic only?') +
  guides(fill = guide_legend(reverse = TRUE)) 

yo = 3

### Southern tonics plots:
level <- 0.95
wt_1_plot <- ggplot(subset(aromanian_tonic, Speaker == "WT_1"), aes(x = F2, y = F1, color = ifelse(VowelType == 2, "/\u0259/", "/\u0268/"))) +
               geom_text(aes(label = VowelInWord), size = yo, alpha = 0.8) +
               scale_y_reverse() + scale_x_reverse() +
               #xlim(14, 11) + ylim(6.5, 3) +
               labs(x = "F2 (Bark)", y = "F1 (Bark)", color = "Vowel") +
               ggtitle(expression(paste("Speaker ", WT[1], sep="")),
                       "Tonic realizations of /\u0259/ and /\u0268/") +
               scale_colour_manual(values=cbPalette) +
               stat_ellipse(type="norm", linetype=1, size=1, alpha = 0.6, level = level) +
               theme_classic() +
               theme(plot.title = element_text(hjust = 0.5),
                     plot.subtitle = element_text(hjust = 0.5)) +
               theme(legend.position = "none")
wt_2_plot <- ggplot(subset(aromanian_tonic, Speaker == "WT_2"), aes(x = F2, y = F1, color = ifelse(VowelType == 2, "/\u0259/", "/\u0268/"))) +
               geom_text(aes(label = VowelInWord), size = yo, alpha = 0.8) +
               scale_y_reverse() + scale_x_reverse() +
               #xlim(14.5, 9) + ylim(7.5, 3.5) +
               labs(x = "F2 (Bark)", y = "F1 (Bark)", color = "Vowel") +
               ggtitle(expression(paste("Speaker ", WT[2], sep="")),
                       "Tonic realizations of /\u0259/ and /\u0268/") +
               scale_colour_manual(values=cbPalette) +
               stat_ellipse(type="norm", linetype=1, size=1, alpha = 0.6, level = level) +
               theme_classic() +
               theme(plot.title = element_text(hjust = 0.5),
                     plot.subtitle = element_text(hjust = 0.5)) +
               theme(legend.position = "none")
mt_1_plot <- ggplot(subset(aromanian_tonic, Speaker == "MT_1"), aes(x = F2, y = F1, color = ifelse(VowelType == 2, "/\u0259/", "/\u0268/"))) +
               geom_text(aes(label = VowelInWord), size = yo, alpha = 0.8) +
               scale_y_reverse() + scale_x_reverse() +
               #xlim(13, 8) + ylim(5.5, 2) +
               labs(x = "F2 (Bark)", y = "F1 (Bark)", color = "Vowel") +
               ggtitle(expression(paste("Speaker ", MT[1], sep="")),
                       "Tonic realizations of /\u0259/ and /\u0268/") +
               scale_colour_manual(values=cbPalette) +
               stat_ellipse(type="norm", linetype=1, size=1, alpha = 0.6, level = level) +
               theme_classic() +
               theme(plot.title = element_text(hjust = 0.5),
                     plot.subtitle = element_text(hjust = 0.5)) +
               theme(legend.position = c(.8, .8),
                     legend.justification = c("center", "center"),
                     legend.box.just = "right",
                     legend.margin = margin(12, 18, 12, 18),
                     legend.background = element_rect(fill=alpha('gray90', 1)))
mt_2_plot <- ggplot(subset(aromanian_tonic, Speaker == "MT_2"), aes(x = F2, y = F1, color = ifelse(VowelType == 2, "/\u0259/", "/\u0268/"))) +
               geom_text(aes(label = VowelInWord), size = yo, alpha = 0.8) +
               scale_y_reverse() + scale_x_reverse() +
               #xlim(13, 9) + ylim(5.0, 3.5) +
               labs(x = "F2 (Bark)", y = "F1 (Bark)", color = "Vowel") +
               ggtitle(expression(paste("Speaker ", MT[2], sep="")),
                       "Tonic realizations of /\u0259/ and /\u0268/") +
               scale_colour_manual(values=cbPalette) +
               stat_ellipse(type="norm", linetype=1, size=1, alpha = 0.6, level = level) +
               theme_classic() +
               theme(plot.title = element_text(hjust = 0.5),
                     plot.subtitle = element_text(hjust = 0.5)) +
               theme(legend.position = "none")

grid.arrange(rbind(cbind(ggplotGrob(wt_1_plot), ggplotGrob(wt_2_plot), size="last"),
                   cbind(ggplotGrob(mt_1_plot), ggplotGrob(mt_2_plot), size="last")))

### Northern tonics plots:
wa_1_plot <- ggplot(subset(aromanian_tonic, Speaker == "WA_1"), aes(x = F2, y = F1, color = ifelse(VowelType == 2, "/\u0259/", "/\u0268/"))) +
               geom_text(aes(label = VowelInWord), size = yo, alpha = 0.8) +
               scale_y_reverse() + scale_x_reverse() +
               #xlim(14.5, 11.5) + ylim(7.5, 3) +
               labs(x = "F2 (Bark)", y = "F1 (Bark)", color = "Vowel") +
               ggtitle(expression(paste("Speaker ", WA[1], sep="")),
                       "Tonic realizations of /\u0259/ and /\u0268/") +
               scale_colour_manual(values=cbPalette) +
               stat_ellipse(type="norm", linetype=1, size=1, alpha = 0.6, level = level) +
               theme_classic() +
               theme(plot.title = element_text(hjust = 0.5),
                     plot.subtitle = element_text(hjust = 0.5)) +
               theme(legend.position = "none")
wa_2_plot <- ggplot(subset(aromanian_tonic, Speaker == "WA_2"), aes(x = F2, y = F1, color = ifelse(VowelType == 2, "/\u0259/", "/\u0268/"))) +
               geom_text(aes(label = VowelInWord), size = yo, alpha = 0.8) +
               scale_y_reverse() + scale_x_reverse() +
               #xlim(15, 12) + ylim(6.5, 4) +
               labs(x = "F2 (Bark)", y = "F1 (Bark)", color = "Vowel") +
               ggtitle(expression(paste("Speaker ", WA[2], sep="")),
                       "Tonic realizations of /\u0259/ and /\u0268/") +
               scale_colour_manual(values=cbPalette) +
               stat_ellipse(type="norm", linetype=1, size=1, alpha = 0.6, level = level) +
               theme_classic() +
               theme(plot.title = element_text(hjust = 0.5),
                     plot.subtitle = element_text(hjust = 0.5)) +
               theme(legend.position = "none")
ma_1_plot <- ggplot(subset(aromanian_tonic, Speaker == "MA_1"), aes(x = F2, y = F1, color = ifelse(VowelType == 2, "/\u0259/", "/\u0268/"))) +
               geom_text(aes(label = VowelInWord), size = yo, alpha = 0.8) +
               scale_y_reverse() + scale_x_reverse() +
               #xlim(14, 6.5) + ylim(6.5, 3.5) +
               labs(x = "F2 (Bark)", y = "F1 (Bark)", color = "Vowel") +
               ggtitle(expression(paste("Speaker ", MA[1], sep="")),
                       "Tonic realizations of /\u0259/ and /\u0268/") +
               scale_colour_manual(values=cbPalette) +
               stat_ellipse(type="norm", linetype=1, size=1, alpha = 0.6, level = level) +
               theme_classic() +
               theme(plot.title = element_text(hjust = 0.5),
                     plot.subtitle = element_text(hjust = 0.5)) +
               theme(legend.position = c(.175, .8),
                     legend.justification = c("center", "center"),
                     legend.box.just = "right",
                     legend.margin = margin(12, 18, 12, 18),
                     legend.background = element_rect(fill=alpha('gray90', 0.75)))
ms_1_plot <- ggplot(subset(aromanian_tonic, Speaker == "MS_1"), aes(x = F2, y = F1, color = ifelse(VowelType == 2, "/\u0259/", "/\u0268/"))) +
               geom_text(aes(label = VowelInWord), size = yo, alpha = 0.8) +
               scale_y_reverse() + scale_x_reverse() +
               #xlim(14, 11.5) + ylim(7, 3) +
               labs(x = "F2 (Bark)", y = "F1 (Bark)", color = "Vowel") +
               ggtitle(expression(paste("Speaker ", MS[1], sep="")),
                       "Tonic realizations of /\u0259/ and /\u0268/") +
               scale_colour_manual(values=cbPalette) +
               stat_ellipse(type="norm", linetype=1, size=1, alpha = 0.6, level = level) +
               theme_classic() +
               theme(plot.title = element_text(hjust = 0.5),
                     plot.subtitle = element_text(hjust = 0.5)) +
               theme(legend.position = "none")

grid.arrange(rbind(cbind(ggplotGrob(wa_1_plot), ggplotGrob(wa_2_plot), size="last"),
                   cbind(ggplotGrob(ma_1_plot), ggplotGrob(ms_1_plot), size="last")))

## tonic mg1

mg_1_plot <- ggplot(subset(aromanian_tonic, Speaker == "MG_1"), aes(x = F2, y = F1, color = ifelse(VowelType == 2, "/\u0259/", "/\u0268/"))) +
               geom_text(aes(label = VowelInWord), size = yo, alpha = 0.8) +
               scale_y_reverse() + scale_x_reverse() +
               #xlim(13.5, 10) + ylim(7.5, 3) +
               labs(x = "F2 (Bark)", y = "F1 (Bark)", color = "Vowel") +
               ggtitle(expression(paste("Speaker ", MG[1], sep="")),
                       "Tonic realizations of /\u0259/ and /\u0268/") +
               scale_colour_manual(values=cbPalette) +
               stat_ellipse(type="norm", linetype=1, size=1, alpha = 0.6, level = level) +
               theme_classic() +
               theme(plot.title = element_text(hjust = 0.5),
                     plot.subtitle = element_text(hjust = 0.5)) +
               theme(legend.position = c(.9, .125),
                     legend.justification = c("center", "center"),
                     legend.box.just = "right",
                     legend.margin = margin(12, 18, 12, 18),
                     legend.background = element_rect(fill=alpha('gray90', 1)))
print(mg_1_plot)

### the mean plots 14:

target_vowels <- c("/\u0259/", "/\u0268/")
other_vowels <- c("/a/", "/e/", "/i/", "/u/", "/o/")

south_w_f1_means <- matrix(vector(), 0, 7)
colnames(south_w_f1_means) <- c(target_vowels, other_vowels)
south_w_f2_means <- matrix(vector(), 0, 7)
colnames(south_w_f2_means) <- c(target_vowels, other_vowels)

south_m_f1_means <- matrix(vector(), 0, 7)
colnames(south_m_f1_means) <- c(target_vowels, other_vowels)
south_m_f2_means <- matrix(vector(), 0, 7)
colnames(south_m_f2_means) <- c(target_vowels, other_vowels)

north_w_f1_means <- matrix(vector(), 0, 7)
colnames(north_w_f1_means) <- c(target_vowels, other_vowels)
north_w_f2_means <- matrix(vector(), 0, 7)
colnames(north_w_f2_means) <- c(target_vowels, other_vowels)

north_m_f1_means <- matrix(vector(), 0, 7)
colnames(north_m_f1_means) <- c(target_vowels, other_vowels)
north_m_f2_means <- matrix(vector(), 0, 7)
colnames(north_m_f2_means) <- c(target_vowels, other_vowels)

for (speaker in unique(aromanian$Speaker)) {
  current_subset <- subset(aromanian, Speaker == speaker)
  other_subset <- subset(all_vowels, Speaker == speaker)
  f1_mean <- mean(current_subset$F1)
  f2_mean <- mean(current_subset$F2)
  region_and_sex <- c(as.vector(unique(current_subset$Region)),
                      as.vector(unique(current_subset$Sex)))
  if (all(region_and_sex == c("South", "Female"))) {
    new_row <- c()
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
                 mean(subset(current_subset, VowelType == 2)$F2),
                 mean(subset(current_subset, VowelType == 1)$F2))
    for (vowel in other_vowels) {
      new_row <- c(new_row, mean(subset(other_subset, VowelType == substr(vowel, 2, 2))$F2))
    }
    south_w_f2_means <- rbind(south_w_f2_means, new_row)
  } else if (all(region_and_sex == c("South", "Male"))) {
    new_row <- c()
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
                 mean(subset(current_subset, VowelType == 2)$F2),
                 mean(subset(current_subset, VowelType == 1)$F2))
    for (vowel in other_vowels) {
      new_row <- c(new_row, mean(subset(other_subset, VowelType == substr(vowel, 2, 2))$F2))
    }
    south_m_f2_means <- rbind(south_m_f2_means, new_row)
  } else if (all(region_and_sex == c("North", "Female"))) {
    new_row <- c()
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
                 mean(subset(current_subset, VowelType == 2)$F2),
                 mean(subset(current_subset, VowelType == 1)$F2))
    for (vowel in other_vowels) {
      new_row <- c(new_row, mean(subset(other_subset, VowelType == substr(vowel, 2, 2))$F2))
    }
    north_w_f2_means <- rbind(north_w_f2_means, new_row)
  } else if (all(region_and_sex == c("North", "Male"))) {
    new_row <- c()
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

for (vowel in c(target_vowels, other_vowels)) {
  south_info <- rbind(south_info, c(vowel, "Female", mean(south_w_f1_means[, vowel]), mean(south_w_f2_means[, vowel])))
  south_info <- rbind(south_info, c(vowel, "Male", mean(south_m_f1_means[, vowel]), mean(south_m_f2_means[, vowel])))
}

for (vowel in c(target_vowels, other_vowels)) {
  north_info <- rbind(north_info, c(vowel, "Female", mean(north_w_f1_means[, vowel]), mean(north_w_f2_means[, vowel])))
  north_info <- rbind(north_info, c(vowel, "Male", mean(north_m_f1_means[, vowel]), mean(north_m_f2_means[, vowel])))
}

south_info <- data.frame(south_info)
north_info <- data.frame(north_info)
south_info <- mutate(south_info, Vowel = as.factor(Vowel), Sex = as.factor(Sex), F1 = as.numeric(F1), F2 = as.numeric(F2))
north_info <- mutate(north_info, Vowel = as.factor(Vowel), Sex = as.factor(Sex), F1 = as.numeric(F1), F2 = as.numeric(F2))

south_w_plot <- ggplot(subset(south_info, Sex == "Female"), aes(x = F2, y = F1)) +
                  geom_text(aes(label = Vowel), size = 3.5) +
                  scale_x_reverse() + scale_y_reverse() +
                  labs(x = "F2 (Bark)", y = "F1 (Bark)") +
                  ggtitle(expression(paste(WT[1], ", ", WT[2], ", ", WK[1], ", ", WK[2], ", ", WK[3], sep="")),
                          "Average Southern-female vowel space") +
                  theme_classic() +
                  theme(plot.title = element_text(hjust = 0.5),
                        plot.subtitle = element_text(hjust = 0.5))

south_m_plot <- ggplot(subset(south_info, Sex == "Male"), aes(x = F2, y = F1)) +
                  geom_text(aes(label = Vowel), size = 3.5) +
                  scale_x_reverse() + scale_y_reverse() +
                  labs(x = "F2 (Bark)", y = "F1 (Bark)") +
                  ggtitle(expression(paste(MT[1], ", ", MT[2], sep="")),
                          "Average Southern-male vowel space") +
                  theme_classic() +
                  theme(plot.title = element_text(hjust = 0.5),
                        plot.subtitle = element_text(hjust = 0.5))

north_w_plot <- ggplot(subset(north_info, Sex == "Female"), aes(x = F2, y = F1)) +
                  geom_text(aes(label = Vowel), size = 3.5) +
                  scale_x_reverse() + scale_y_reverse() +
                  labs(x = "F2 (Bark)", y = "F1 (Bark)") +
                  ggtitle(expression(paste(WA[1], ", ", WA[2], ", ", WA[3], sep="")),
                          "Average Northern-female vowel space") +
                  theme_classic() +
                  theme(plot.title = element_text(hjust = 0.5),
                        plot.subtitle = element_text(hjust = 0.5))

north_m_plot <- ggplot(subset(north_info, Sex == "Male"), aes(x = F2, y = F1)) +
                  geom_text(aes(label = Vowel), size = 3.5) +
                  scale_x_reverse() + scale_y_reverse() +
                  labs(x = "F2 (Bark)", y = "F1 (Bark)") +
                  ggtitle(expression(paste(MA[1], ", ", MA[2], ", ", MG[1], ", ", MS[1], sep="")),
                          "Average Northern-male vowel space") +
                  theme_classic() +
                  theme(plot.title = element_text(hjust = 0.5),
                        plot.subtitle = element_text(hjust = 0.5))

grid.arrange(rbind(cbind(ggplotGrob(south_w_plot), ggplotGrob(north_w_plot), size="last"),
                   cbind(ggplotGrob(south_m_plot), ggplotGrob(north_m_plot), size="last")))

### the mean plots tonic 9:

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
  region_and_sex <- c(as.vector(unique(current_subset$Region)),
                      as.vector(unique(current_subset$Sex)))
  if (all(region_and_sex == c("South", "Female"))) {
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
  } else if (all(region_and_sex == c("South", "Male"))) {
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
  } else if (all(region_and_sex == c("North", "Female"))) {
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
  } else if (all(region_and_sex == c("North", "Male"))) {
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

customPalette <- c("black", "#009E73")

south_w_plot <- ggplot(subset(south_info, Sex == "Female"), aes(x = F2, y = F1, color = ifelse(!is.na(str_match(Vowel, "tonic")), "yes", "no"))) +
  geom_text(aes(label = substring(Vowel, 1, 3)), size = 3.5) +
  scale_x_reverse() + scale_y_reverse() +
  labs(x = "F2 (Bark)", y = "F1 (Bark)") +
  ggtitle(expression(paste(WT[1], ", ", WT[2], sep="")),
          "Average Southern-female vowel space") +
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
          "Average Southern-male vowel space") +
  scale_color_manual(values=customPalette) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  theme(legend.position = "none")

north_w_plot <- ggplot(subset(north_info, Sex == "Female"), aes(x = F2, y = F1, color = ifelse(!is.na(str_match(Vowel, "tonic")), "yes", "no"))) +
  geom_text(aes(label = substring(Vowel, 1, 3)), size = 3.5) +
  scale_x_reverse() + scale_y_reverse() +
  labs(x = "F2 (Bark)", y = "F1 (Bark)", color = "Tonic only?") +
  ggtitle(expression(paste(WA[1], ", ", WA[2], sep="")),
          "Average Northern-female vowel space") +
  scale_color_manual(values=customPalette) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  theme(legend.position = c(.75, .3),
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
          "Average Northern-male vowel space") +
  scale_color_manual(values=customPalette) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  theme(legend.position = "none")

grid.arrange(rbind(cbind(ggplotGrob(south_w_plot), ggplotGrob(north_w_plot), size="last"),
                   cbind(ggplotGrob(south_m_plot), ggplotGrob(north_m_plot), size="last")))


## weiter gehts
aromanian_tonic_fresh <- subset(aromanian, Stress == "yes")

first_cont <- table(aromanian$VowelType, aromanian$Stress)

second_cont <- table(aromanian_tonic_fresh$VowelType, aromanian_tonic_fresh$POS)

third_cont <- table(aromanian_tonic_fresh$VowelType, aromanian_tonic_fresh$LeftCoArt)

fourth_cont <- table(aromanian_tonic_fresh$VowelType, aromanian_tonic_fresh$RightCoArt)

### POS and CoArt beispiele

convert_context <- function(context) {
  if (context %in% c("^_", "_$")) {
    return("C1")
  } else if (context %in% c("a_", "_a")) {
    return("C2")
  } else if (context %in% c("e_", "_e")) {
    return("C3")
  } else if (context %in% c("i_", "_i")) {
    return("C4")
  } else if (context %in% c("o_", "_o")) {
    return("C5")
  } else if (context %in% c("u_", "_u")) {
    return("C6")
  } else if (context %in% c("@_", "_@")) {
    return("C7")
  } else if (context %in% c("m_", "_m")) {
    return("C8")
  } else if (context %in% c("n_", "_n")) {
    return("C9")
  } else if (context %in% c("N_", "_N")) {
    return("C10")
  } else if (context %in% c("J_", "_J")) {
    return("C11")
  } else if (context %in% c("l_", "_l")) {
    return("C12")
  } else if (context %in% c("L_", "_L")) {
    return("C13")
  } else if (context %in% c("r_", "_r")) {
    return("C14")
  } else if (context %in% c("h_", "_h")) {
    return("C15")
  } else if (context %in% c("b/p_", "_b/p")) {
    return("C16")
  } else if (context %in% c("v/f_", "_v/f")) {
    return("C17")
  } else if (context %in% c("D/T_", "_D/T")) {
    return("C18")
  } else if (context %in% c("G/R/x_", "_G/R/x")) {
    return("C19")
  } else if (context == "j/c/C_") {
    return("P1")
  } else if (context == "g/k_") {
    return("P2")
  } else if (context == "d/t_") {
    return("P3")
  } else if (context == "z/s/dz/ts_") {
    return("P4")
  } else if (context == "Z/S/dZ/tS_") {
    return("P5")
  } else if (context == "_g/c/k") {
    return("S1")
  } else if (context == "_j/C") {
    return("S2")
  } else if (context == "_z/s") {
    return("S3")
  } else if (context == "_Z/S") {
    return("S4")
  } else if (context == "_d/t/dz/ts/dZ/tS") {
    return("S5")
  } else {
    print(context)
    return(NULL)
  }
}

wt_1_plot_POS <- ggplot(subset(aromanian_tonic, Speaker == "WT_1"), aes(x = F2, y = F1, color = ifelse(VowelType == 2, "/\u0259/", "/\u0268/"))) +
  geom_text(aes(label = POS), size = yo, alpha = 0.9) +
  #scale_y_reverse() + scale_x_reverse() +
  xlim(14, 11) + ylim(6.5, 3) +
  labs(x = "F2 (Bark)", y = "F1 (Bark)", color = "Vowel") +
  ggtitle(expression(paste("Speaker ", WT[1], sep="")),
          "Tonic realizations of /\u0259/ and /\u0268/") +
  scale_colour_manual(values=cbPalette) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  theme(legend.position = c(0.975, 0.8),
        legend.justification = c("center", "center"),
        legend.box.just = "right",
        legend.margin = margin(6, 9, 6, 9),
        legend.background = element_rect(fill=alpha('gray90', 1)))
#print(wt_1_plot_POS)

wt_1_plot_CA <- ggplot(subset(aromanian_tonic, Speaker == "WT_1"), aes(x = F2, y = F1, color = ifelse(VowelType == 2, "/\u0259/", "/\u0268/"))) +
  geom_text(aes(label = paste0(sapply(LeftCoArt, convert_context), "-", sapply(RightCoArt, convert_context))), size = yo, alpha = 0.9) +
  #scale_y_reverse() + scale_x_reverse() +
  xlim(14, 11) + ylim(6.5, 3) +
  labs(x = "F2 (Bark)", y = "F1 (Bark)", color = "Vowel") +
  ggtitle(expression(paste("Speaker ", WT[1], sep="")),
          "Tonic realizations of /\u0259/ and /\u0268/") +
  scale_colour_manual(values=cbPalette) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  theme(legend.position = "none")
#print(wt_1_plot_CA)

grid.arrange(cbind(ggplotGrob(wt_1_plot_POS), ggplotGrob(wt_1_plot_CA), size="last"))

## P4V stuff

aromanian_tonic_wop4v <- subset(aromanian_tonic, !(LeftCoArt == "z/s/dz/ts_" & POS == "V"))
complementary_here <- subset(aromanian_tonic, LeftCoArt == "z/s/dz/ts_" & POS == "V")
table(aromanian_tonic_wop4v$Speaker, aromanian_tonic_wop4v$VowelType) # this is fine

## wop4v manova
factors <- c("VowelType", "VowelType+D.")
complete_cmb_list <- list()
for (speaker in unique(aromanian_tonic_wop4v$Speaker)) {
  print(speaker)
  current_subset <- subset(aromanian_tonic_wop4v, Speaker == speaker)
  cmb <- list()
  
  cmb[[factors[1]]] <- manova(cbind(F1, F2, F3) ~ VowelType, data = current_subset)
  cmb[[factors[2]]] <- manova(cbind(F1, F2, F3) ~ VowelType + Duration, data = current_subset)
  
  complete_cmb_list[[speaker]] <- cmb
}

wop4v_tonic_pillai_results <- c()
wop4v_tonic_p_results <- c()
for (speaker in unique(aromanian_tonic_wop4v$Speaker)) {
  print(speaker)
  for (fct in factors) {
    summary_mnv <- summary(complete_cmb_list[[speaker]][[fct]])
    wop4v_tonic_pillai <- summary_mnv$stats[, "Pillai"][str_replace_all(fct, "\\+D\\.", "")]
    wop4v_tonic_p <- summary_mnv$stats[, "Pr(>F)"][str_replace_all(fct, "\\+D\\.", "")]
    wop4v_tonic_pillai_results <- c(wop4v_tonic_pillai_results, wop4v_tonic_pillai)
    wop4v_tonic_p_results <- c(wop4v_tonic_p_results, wop4v_tonic_p)
  }
}

wop4v_tonic_pillai_matrix <- matrix(wop4v_tonic_pillai_results, nrow = length(factors), ncol = length(unique(aromanian_tonic_wop4v$Speaker)))
wop4v_tonic_p_matrix <- matrix(wop4v_tonic_p_results, nrow = length(factors), ncol = length(unique(aromanian_tonic_wop4v$Speaker)))

colnames(wop4v_tonic_pillai_matrix) <- unique(aromanian_tonic_wop4v$Speaker)
colnames(wop4v_tonic_p_matrix) <- unique(aromanian_tonic_wop4v$Speaker)

wop4v_tonic_pillai_rownames <- c()
wop4v_tonic_p_rownames <- c() 
for (fct in factors) {
  wop4v_tonic_pillai_rownames <- c(wop4v_tonic_pillai_rownames, paste(fct, "Pillai"))
  wop4v_tonic_p_rownames <- c(wop4v_tonic_p_rownames, paste(fct, "p-value"))
}
rownames(wop4v_tonic_pillai_matrix) <- wop4v_tonic_pillai_rownames
rownames(wop4v_tonic_p_matrix) <- wop4v_tonic_p_rownames

wop4v_tonic_significance_matrix <- wop4v_tonic_p_matrix < 0.05

wop4v_tonic_pillai_matrix_signif <- matrix(nrow = nrow(wop4v_tonic_p_matrix), ncol = ncol(wop4v_tonic_p_matrix))
wop4v_tonic_p_matrix_signif <- matrix(nrow = nrow(wop4v_tonic_p_matrix), ncol = ncol(wop4v_tonic_p_matrix))

colnames(wop4v_tonic_pillai_matrix_signif) <- unique(aromanian_tonic_wop4v$Speaker)
colnames(wop4v_tonic_p_matrix_signif) <- unique(aromanian_tonic_wop4v$Speaker)
rownames(wop4v_tonic_pillai_matrix_signif) <- wop4v_tonic_pillai_rownames
rownames(wop4v_tonic_p_matrix_signif) <- wop4v_tonic_p_rownames

for (i in 1:nrow(wop4v_tonic_p_matrix)) {
  for (j in 1:ncol(wop4v_tonic_p_matrix)) {
    if (!is.na(wop4v_tonic_significance_matrix[i, j]) & wop4v_tonic_significance_matrix[i, j]) {
      wop4v_tonic_pillai_matrix_signif[i, j] = wop4v_tonic_pillai_matrix[i, j]
      wop4v_tonic_p_matrix_signif[i, j] = wop4v_tonic_p_matrix[i, j]
    }
  }
}

write.csv(data.frame(wop4v_tonic_pillai_matrix), "wop4v_tonic_pillai_matrix.csv")

write.csv(data.frame(wop4v_tonic_p_matrix), "wop4v_tonic_p_matrix.csv")

write.csv(data.frame(wop4v_tonic_pillai_matrix_signif), "wop4v_tonic_pillai_matrix_signif.csv")

###

#copy-pasted geom_bar stuff die zweite:

wop4v_tonic_pillai_df <- data.frame(t(wop4v_tonic_pillai_matrix))
colnames(wop4v_tonic_pillai_df) <- factors
wop4v_tonic_pillai_df$Speaker <- factor(rownames(wop4v_tonic_pillai_df), levels = rev(rownames(wop4v_tonic_pillai_df)))
rownames(wop4v_tonic_pillai_df) <- NULL
wop4v_tonic_pillai_likethis_df <- subset(wop4v_tonic_pillai_df, select = c(factors[1], "Speaker"))
gentonic_pillai_likethis_df <- subset(tonic_pillai_likethis_df, select = c("VowelType", "Speaker"))

wop4v_tonic_pillai_likethis_df$WOP4V <- ". without P4 verbs"
gentonic_pillai_likethis_df$WOP4V <- ". all"
wop4v_tonic_pillai_df <- rbind(gentonic_pillai_likethis_df, wop4v_tonic_pillai_likethis_df)
wop4v_tonic_pillai_df$WOP4V <- factor(wop4v_tonic_pillai_df$WOP4V, levels = c(". without P4 verbs", ". all"))

wop4v_tonic_speaker_exps <- expression(WT[1], WT[2], MT[1], MT[2],
                                       WA[1], WA[2], MA[1], MG[1], MS[1])

#cbPalette <- c("#56B4E9", "#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00", "chartreuse3", "#999999")
customPalette <- c("chartreuse3", "#009E73")

ggplot(wop4v_tonic_pillai_df, aes(x = Speaker, y = VowelType, fill = WOP4V)) +
  geom_bar(position="dodge", stat="identity", width = 0.8) +
  labs(x = "Speaker", y = "Pillai score between /\u0259/ and /\u0268/") +
  ylim(0, 1) +
  coord_flip() +
  scale_x_discrete(breaks = unique(aromanian$Speaker),
                   labels = speaker_exps) +
  scale_fill_manual(values=customPalette) +
  theme_classic() +
  theme(legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(12, 18, 12, 18),
        legend.background = element_rect(fill=alpha('gray90', 1))) +
  labs(fill = 'Tonic only .') +
  guides(fill = guide_legend(reverse = TRUE)) 

## can do for the four manovas:

customPalette <- c("darkseagreen3", "#009E73")
mg1_stress_plot <- ggplot(subset(aromanian, Speaker == "MG_1"), aes(x = F2, y = F1, color = ifelse(Stress == "yes", "tonic", "non-tonic"))) +
                     geom_point(alpha = 0.8) +
                     scale_x_reverse() + scale_y_reverse() +
                     stat_ellipse(type="norm", linetype=1, size=1, alpha = 0.6, level = level) +
                     theme_classic() +
                     scale_color_manual(values=customPalette) +
                     labs(x = "F2 (Bark)", y = "F1 (Bark)", color = "Stress") +
                     theme(legend.position = "none") +
                     guides(color = guide_legend(reverse = TRUE)) +
                     ggtitle(expression(paste("Northern speaker ", MG[1], sep="")),
                             "Realizations of /\u0259/ and /\u0268/") +
                     theme(plot.title = element_text(hjust = 0.5),
                           plot.subtitle = element_text(hjust = 0.5))
ms1_stress_plot <- ggplot(subset(aromanian, Speaker == "MS_1"), aes(x = F2, y = F1, color = ifelse(Stress == "yes", "tonic", "non-tonic"))) +
                     geom_point(alpha = 0.8) +
                     scale_x_reverse() + scale_y_reverse() +
                     stat_ellipse(type="norm", linetype=1, size=1, alpha = 0.6, level = level) +
                     theme_classic() +
                     scale_color_manual(values=customPalette) +
                     labs(x = "F2 (Bark)", y = "F1 (Bark)", color = "Stress") +
                     theme(legend.position = "none") +
                     guides(color = guide_legend(reverse = TRUE)) +
                     ggtitle(expression(paste("Northern speaker ", MS[1], sep="")),
                             "Realizations of /\u0259/ and /\u0268/") +
                     theme(plot.title = element_text(hjust = 0.5),
                           plot.subtitle = element_text(hjust = 0.5))
wa2_stress_plot <- ggplot(subset(aromanian, Speaker == "WA_2"), aes(x = F2, y = F1, color = ifelse(Stress == "yes", "tonic", "non-tonic"))) +
                     geom_point(alpha = 0.8) +
                     scale_x_reverse() + scale_y_reverse() +
                     stat_ellipse(type="norm", linetype=1, size=1, alpha = 0.6, level = level) +
                     theme_classic() +
                     scale_color_manual(values=customPalette) +
                     labs(x = "F2 (Bark)", y = "F1 (Bark)", color = "Stress") +
                     theme(legend.position = "none") +
                     guides(color = guide_legend(reverse = TRUE)) +
                     ggtitle(expression(paste("Northern speaker ", WA[2], sep="")),
                             "Realizations of /\u0259/ and /\u0268/") +
                     theme(plot.title = element_text(hjust = 0.5),
                           plot.subtitle = element_text(hjust = 0.5)) +
                     theme(legend.position = c(1.1, 1.0),
                           legend.justification = c("right", "top"),
                           legend.box.just = "right",
                           legend.margin = margin(4, 6, 4, 6),
                           legend.background = element_rect(fill=alpha('gray90', 1)))
mt1_stress_plot <- ggplot(subset(aromanian, Speaker == "MT_1"), aes(x = F2, y = F1, color = ifelse(Stress == "yes", "tonic", "non-tonic"))) +
                     geom_point(alpha = 0.8) +
                     scale_x_reverse() + scale_y_reverse() +
                     stat_ellipse(type="norm", linetype=1, size=1, alpha = 0.6, level = level) +
                     theme_classic() +
                     scale_color_manual(values=customPalette) +
                     labs(x = "F2 (Bark)", y = "F1 (Bark)", color = "Stress") +
                     theme(legend.position = "none") +
                     guides(color = guide_legend(reverse = TRUE)) +
                     ggtitle(expression(paste("Southern speaker ", MT[1], sep="")),
                             "Realizations of /\u0259/ and /\u0268/") +
                     theme(plot.title = element_text(hjust = 0.5),
                           plot.subtitle = element_text(hjust = 0.5))

grid.arrange(rbind(cbind(ggplotGrob(mg1_stress_plot), ggplotGrob(ms1_stress_plot), size="last"),
                   cbind(ggplotGrob(wa2_stress_plot), ggplotGrob(mt1_stress_plot), size="last")))


###
customPalette <- c("#D55E00", "chartreuse3")

mg1_prost_plot <- ggplot(subset(aromanian, Speaker == "MG_1" & Stress %in% c("pre", "post")), aes(x = F2, y = F1, color = Stress)) +
                    geom_point(alpha = 0.8) +
                    scale_x_reverse() + scale_y_reverse() +
                    stat_ellipse(type="norm", linetype=1, size=1, alpha = 0.6) +
                    theme_classic() +
                    scale_color_manual(values=customPalette) +
                    labs(x = "F2 (Bark)", y = "F1 (Bark)", color = "Stress") +
                    theme(legend.position = "none") +
                    guides(color = guide_legend(reverse = TRUE)) +
                    ggtitle(expression(paste("Northern speaker ", MG[1], sep="")),
                            "Realizations of /\u0259/ and /\u0268/") +
                    theme(plot.title = element_text(hjust = 0.5),
                          plot.subtitle = element_text(hjust = 0.5))
ms1_prost_plot <- ggplot(subset(aromanian, Speaker == "MS_1" & Stress %in% c("pre", "post")), aes(x = F2, y = F1, color = Stress)) +
                    geom_point(alpha = 0.8) +
                    scale_x_reverse() + scale_y_reverse() +
                    stat_ellipse(type="norm", linetype=1, size=1, alpha = 0.6) +
                    theme_classic() +
                    scale_color_manual(values=customPalette) +
                    labs(x = "F2 (Bark)", y = "F1 (Bark)", color = "Stress") +
                    theme(legend.position = "none") +
                    guides(color = guide_legend(reverse = TRUE)) +
                    ggtitle(expression(paste("Northern speaker ", MS[1], sep="")),
                            "Realizations of /\u0259/ and /\u0268/") +
                    theme(plot.title = element_text(hjust = 0.5),
                          plot.subtitle = element_text(hjust = 0.5))
wa2_prost_plot <- ggplot(subset(aromanian, Speaker == "WA_1" & Stress %in% c("pre", "post")), aes(x = F2, y = F1, color = Stress)) +
                    geom_point(alpha = 0.8) +
                    scale_x_reverse() + scale_y_reverse() +
                    stat_ellipse(type="norm", linetype=1, size=1, alpha = 0.6) +
                    theme_classic() +
                    scale_color_manual(values=customPalette) +
                    labs(x = "F2 (Bark)", y = "F1 (Bark)", color = "Stress") +
                    theme(legend.position = "none") +
                    guides(color = guide_legend(reverse = TRUE)) +
                    ggtitle(expression(paste("Northern speaker ", WA[2], sep="")),
                            "Realizations of /\u0259/ and /\u0268/") +
                    theme(plot.title = element_text(hjust = 0.5),
                          plot.subtitle = element_text(hjust = 0.5)) +
                    theme(legend.position = c(1.1, 1.0),
                          legend.justification = c("right", "top"),
                          legend.box.just = "right",
                          legend.margin = margin(4, 6, 4, 6),
                          legend.background = element_rect(fill=alpha('gray90', 1)))
mt1_prost_plot <- ggplot(subset(aromanian, Speaker == "MT_1" & Stress %in% c("pre", "post")), aes(x = F2, y = F1, color = Stress)) +
                    geom_point(alpha = 0.8) +
                    scale_x_reverse() + scale_y_reverse() +
                    stat_ellipse(type="norm", linetype=1, size=1, alpha = 0.6) +
                    theme_classic() +
                    scale_color_manual(values=customPalette) +
                    labs(x = "F2 (Bark)", y = "F1 (Bark)", color = "Stress") +
                    theme(legend.position = "none") +
                    guides(color = guide_legend(reverse = TRUE)) +
                    ggtitle(expression(paste("Southern speaker ", MT[1], sep="")),
                            "Realizations of /\u0259/ and /\u0268/") +
                    theme(plot.title = element_text(hjust = 0.5),
                          plot.subtitle = element_text(hjust = 0.5))

grid.arrange(rbind(cbind(ggplotGrob(mg1_prost_plot), ggplotGrob(ms1_prost_plot), size="last"),
                   cbind(ggplotGrob(wa2_prost_plot), ggplotGrob(mt1_prost_plot), size="last")))

# okblabla
'
ggplot(subset(aromanian, Speaker == "MS_1" & RightCoArt %in% c("_$", "_b/p", "_d/t/dz/ts/dZ/tS", "_g/c/k", "_l", "_m", "_n", "_z/s")), aes(x = F2, y = F1, color = RightCoArt)) +
  geom_point(alpha = 0.6) +
  stat_ellipse(type="norm", linetype=1, size=1, alpha = 0.6, level = 0.01) +
  scale_x_reverse() + scale_y_reverse()
'

mg1_sub <- subset(aromanian, Speaker == "MG_1")
ms1_sub <- subset(aromanian, Speaker == "MS_1")
wa2_sub <- subset(aromanian, Speaker == "WA_2")
mt1_sub <- subset(aromanian, Speaker == "MT_1")

### sifting MG1 for 0 < 3

# voweltype, stress, pos OK
# leftcoart: a 1, jcC 1, o 1

mg1_sub <- subset(mg1_sub, !(LeftCoArt %in% c("a_", "j/c/C_", "o_")))

# rightcoart: DT 2, e 1

mg1_sub <- subset(mg1_sub, !(RightCoArt %in% c("_D/T", "_e")))

# pos: m 2

mg1_sub <- subset(mg1_sub, !(POS %in% c("M")))

# leftcoart: i 2

mg1_sub <- subset(mg1_sub, !(LeftCoArt %in% c("i_")))

# all OK

#####

### sifting MS1 for 0 < 3

# voweltype, stress OK
# pos: a 2, m 2

ms1_sub <- subset(ms1_sub, !(POS %in% c("A", "M")))

# leftcoart: @ 1, a 2, l 2

ms1_sub <- subset(ms1_sub, !(LeftCoArt %in% c("@_", "a_", "l_")))

# rightcoart: @ 1, e 1, GRx 2, jC 2, u 2, ZS 1

ms1_sub <- subset(ms1_sub, !(RightCoArt %in% c("_@", "_e", "_G/R/x", "_j/C", "_u", "_Z/S")))

# all OK

#####

### sifting WA2 for 0 < 3

# voweltype, stress OK
# pos: a 1

wa2_sub <- subset(wa2_sub, !(POS %in% c("A")))

# leftcoart: ^ 1

wa2_sub <- subset(wa2_sub, !(LeftCoArt %in% c("^_")))

# rightcoart: a 2, e 2, i 1, jC 1, L 1, N 2, o 1

wa2_sub <- subset(wa2_sub, !(RightCoArt %in% c("_a", "_e", "_i", "_j/C", "_L", "_N", "_o")))

# leftcoart: ZSdZtS 2

wa2_sub <- subset(wa2_sub, !(LeftCoArt %in% c("Z/S/dZ/tS_")))

# pos: r 2, p 2

wa2_sub <- subset(wa2_sub, !(POS %in% c("R", "P")))

# leftcoart: l 2

wa2_sub <- subset(wa2_sub, !(LeftCoArt %in% c("l_")))

# all OK

#####

### sifting MT1 for 0 < 3

# voweltype, stress, pos OK
# leftcoart: ^ 2, a 1, i 2

mt1_sub <- subset(mt1_sub, !(LeftCoArt %in% c("^_", "a_", "i_")))

# rightcoart: GRx 1, u 1

mt1_sub <- subset(mt1_sub, !(RightCoArt %in% c("_G/R/x", "_u")))

# all OK

#####

mg1_mod <- lm(cbind(F1, F2) ~ VowelType + POS + LeftCoArt + RightCoArt + Stress, data = mg1_sub)
ms1_mod <- lm(cbind(F1, F2) ~ VowelType + POS + LeftCoArt + RightCoArt + Stress, data = ms1_sub)
wa2_mod <- lm(cbind(F1, F2) ~ VowelType + POS + LeftCoArt + RightCoArt + Stress, data = wa2_sub)
mt1_mod <- lm(cbind(F1, F2) ~ VowelType + POS + LeftCoArt + RightCoArt + Stress, data = mt1_sub)

##########

###yoyodaemmerens.

ma_1_plot_f3_before <- ggplot(subset(aromanian_tonic, Speaker == "MA_1"), aes(x = F2, y = F3, color = ifelse(VowelType == 2, "/\u0259/", "/\u0268/"))) +
  geom_text(aes(label = VowelInWord), size = yo, alpha = 0.8) +
  #scale_y_reverse() + scale_x_reverse() +
  xlim(14.5, 6) + ylim(16, 13) +
  labs(x = "F2 (Bark)", y = "F3 (Bark)", color = "Vowel") +
  ggtitle(expression(paste("Speaker ", MA[1], sep="")),
          "All tonic realizations of /\u0259/ and /\u0268/") +
  scale_colour_manual(values=cbPalette) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  theme(legend.position = "none")
ma_1_plot_f3_after <- ggplot(subset(aromanian_tonic_wop4v, Speaker == "MA_1"), aes(x = F2, y = F3, color = ifelse(VowelType == 2, "/\u0259/", "/\u0268/"))) +
  geom_text(aes(label = VowelInWord), size = yo, alpha = 0.8) +
  #scale_y_reverse() + scale_x_reverse() +
  xlim(14.5, 6) + ylim(16, 13) +
  labs(x = "F2 (Bark)", y = "F3 (Bark)", color = "Vowel") +
  ggtitle(expression(paste("Speaker ", MA[1], sep="")),
          "Non-P4-verb tonic realizations of /\u0259/ and /\u0268/") +
  scale_colour_manual(values=cbPalette) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  theme(legend.position = c(.8, .2),
        legend.justification = c("center", "center"),
        legend.box.just = "right",
        legend.margin = margin(12, 18, 12, 18),
        legend.background = element_rect(fill=alpha('gray90', 1)))
grid.arrange(cbind(ggplotGrob(ma_1_plot_f3_before), ggplotGrob(ma_1_plot_f3_after), size="last"))


##### experimental

customPalette <- c("#D55E00", "chartreuse3")

mg1_exp_plot <- ggplot(subset(aromanian, Speaker == "MG_1"), aes(x = F2, y = F1, color = RightCoArt)) +
                    geom_text(aes(label = VowelInWord), alpha = 0.5) +
                    scale_x_reverse() + scale_y_reverse() +
                    stat_ellipse(type="norm", linetype=1, size=1, level = 0.05, alpha = 1) +
                    theme_classic() +
                    labs(x = "F2 (Bark)", y = "F1 (Bark)", color = "RightCoArt") +
                    guides(color = guide_legend(reverse = TRUE)) +
                    ggtitle(expression(paste("Northern speaker ", MG[1], sep="")),
                            "Realizations of /\u0259/ and /\u0268/") +
                    theme(plot.title = element_text(hjust = 0.5),
                          plot.subtitle = element_text(hjust = 0.5))

print(mg1_exp_plot)

ms1_prost_plot <- ggplot(subset(aromanian, Speaker == "MS_1" & Stress %in% c("pre", "post")), aes(x = F2, y = F1, color = Stress)) +
                    geom_point(alpha = 0.8) +
                    scale_x_reverse() + scale_y_reverse() +
                    stat_ellipse(type="norm", linetype=1, size=1, alpha = 0.6) +
                    theme_classic() +
                    scale_color_manual(values=customPalette) +
                    labs(x = "F2 (Bark)", y = "F1 (Bark)", color = "Stress") +
                    theme(legend.position = "none") +
                    guides(color = guide_legend(reverse = TRUE)) +
                    ggtitle(expression(paste("Northern speaker ", MS[1], sep="")),
                            "Realizations of /\u0259/ and /\u0268/") +
                    theme(plot.title = element_text(hjust = 0.5),
                          plot.subtitle = element_text(hjust = 0.5))
wa2_prost_plot <- ggplot(subset(aromanian, Speaker == "WA_1" & Stress %in% c("pre", "post")), aes(x = F2, y = F1, color = Stress)) +
                    geom_point(alpha = 0.8) +
                    scale_x_reverse() + scale_y_reverse() +
                    stat_ellipse(type="norm", linetype=1, size=1, alpha = 0.6) +
                    theme_classic() +
                    scale_color_manual(values=customPalette) +
                    labs(x = "F2 (Bark)", y = "F1 (Bark)", color = "Stress") +
                    theme(legend.position = "none") +
                    guides(color = guide_legend(reverse = TRUE)) +
                    ggtitle(expression(paste("Northern speaker ", WA[2], sep="")),
                            "Realizations of /\u0259/ and /\u0268/") +
                    theme(plot.title = element_text(hjust = 0.5),
                          plot.subtitle = element_text(hjust = 0.5)) +
                    theme(legend.position = c(1.1, 1.0),
                          legend.justification = c("right", "top"),
                          legend.box.just = "right",
                          legend.margin = margin(4, 6, 4, 6),
                          legend.background = element_rect(fill=alpha('gray90', 1)))
mt1_prost_plot <- ggplot(subset(aromanian, Speaker == "MT_1" & Stress %in% c("pre", "post")), aes(x = F2, y = F1, color = Stress)) +
                    geom_point(alpha = 0.8) +
                    scale_x_reverse() + scale_y_reverse() +
                    stat_ellipse(type="norm", linetype=1, size=1, alpha = 0.6) +
                    theme_classic() +
                    scale_color_manual(values=customPalette) +
                    labs(x = "F2 (Bark)", y = "F1 (Bark)", color = "Stress") +
                    theme(legend.position = "none") +
                    guides(color = guide_legend(reverse = TRUE)) +
                    ggtitle(expression(paste("Southern speaker ", MT[1], sep="")),
                            "Realizations of /\u0259/ and /\u0268/") +
                    theme(plot.title = element_text(hjust = 0.5),
                          plot.subtitle = element_text(hjust = 0.5))
