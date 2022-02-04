#' Simulate outcome of competition between 2 immune escape strategies
library(ggpubr)
library(ggsci)
library(AtmRay)
library(dplyr)
library(cowplot)
library(viridis)
library(colorspace)
library(pals)

theme_set(theme_pubr())
### PARAMETERS ###
ANT1_IDX <- 1
ANT2_IDX <- 2
P1_IDX <- 3
P2_IDX <- 4
S1_IDX <- 5
S2_IDX <- 6
default_r <- 0.2

calc_d <- function(a, p, s, d =0.1){
  #' Calculate death rate
  #'
  #' @param a Antigenicity (0-1)
  #' @param p Protection from blockade (0-1)
  #' @param s Suppression (0-1)
  #' @param d Intrinsic death rate (0-1)
  #' @return total death rate

  if(a == 0){
    kill <- 0
  }else{
    kill <- (1 - s/a)*(1 - p)*a
    if(kill < 0){
      kill <- 0
    }
  }
  return(kill + d )
}

compete <- function(p1_a, p1_p, p1_s, p1_name, p2_a, p2_p, p2_s, p2_name, r=default_r){
  #' Determine outcome of compeition between population p1 and population p2
  #'
  #' @param p1_a Antigenicity of p1 (0-1)
  #' @param p1_p Protection of p1 (0-1)
  #' @param p1_2 Suppression of p1 (0-1)
  #' @param p1_name Name of p1
  #' @param p2_a Antigenicity of p2 (0-1)
  #' @param p2_p Protection of p2 (0-1)
  #' @param p2_2 Suppression of p2 (0-1)
  #' @param p2_name Name of p2
  #' @param default_r Growth rate
  #' @return outcome of competition

  d1 <- calc_d(p1_a, p = p1_p, s = p1_s) ### pop 1
  d2 <- calc_d(p2_a, p = p2_p, s = p2_s) ### pop 2


  p1_net_r <- r - d1
  p2_net_r <- r - d2

  a12 <- (r-d2)/((r-d1)*(1-0.5*p2_s))
  a21 <- (r-d1)/((r-d2)*(1-0.5*p1_s))

  if(p1_net_r <= 0 & p2_net_r <= 0){
    # result <- "Both extinct"
    result <- "X"
  }else if(p1_net_r > 0 & p2_net_r <= 0){
    ### Case 3a
    result <- p1_name
  }else if(p1_net_r <= 0 & p2_net_r > 0){
    ### Case 3b
    result <- p2_name
  }else if(a12 <= 1 & a21 <= 1 ){
    ### Case 1
    # result <- "Coexistence"
    result <- "Co"
  }else if(a12 > 1 & a21 > 1){
    ### Case 2
    # result <- "Initial size"
    result <- "N"
  }else if((a12 <= 1 & a21 >= 1) | (p1_net_r > 0 & p2_net_r < 0)){
    ### Case 3a
    result <- p1_name
  }else if((a12 >= 1 & a21 <= 1) | (p1_net_r < 0 & p2_net_r > 0)){
    ### Case 3b
    result <- p2_name
  }else{
    result <- NA
  }


  reults_df <- data.frame(p1_name, p2_name, a12, a21, p1_net_r, p2_net_r, result)
  reults_df[paste(p1_name, "Ant", sep="_")] <- signif(p1_a, 3)
  reults_df[paste(p1_name, "Pro", sep="_")] <- signif(p1_p, 3)
  reults_df[paste(p1_name, "Sup", sep="_")] <- signif(p1_s, 3)
  reults_df[paste(p2_name, "Ant", sep="_")] <- signif(p2_a, 3)
  reults_df[paste(p2_name, "Pro", sep="_")] <- signif(p2_p, 3)
  reults_df[paste(p2_name, "Sup", sep="_")] <- signif(p2_s, 3)


  return(reults_df)
}


calc_max_size <- function(a, p, s, r=default_r, d =0.1, k = 100){
  #'Determine size at which birth = death. This is the population's largest possible size
  #'
  #' @param a Antigenicity (0-1)
  #' @param p Protection from blockade (0-1)
  #' @param s Suppression (0-1)
  #' @param r Division rate
  #' @param d Intrinsic death rate (0-1)
  #' @param k Carrying capacity
  #' @return N Maximum population size

  death_r <- calc_d(a, p, s, d)
  N <- (1/(0.5*s-1))*((death_r*k)/r - k)

  if(N < 0){
  N <- 0
  }

  return(N)
}
calc_max_size_vec <- Vectorize(calc_max_size)
facet_font_size <- 9

#### Number of grid points in the parameter sweep grid ####
nant <- 20 # Number of antigen levels per strategy
nstrat <- 4 # Number of levels for each strategy.

#### Build parameter grid ####
ant1 <- seq(0, 1, 1/nant) # Antigenicity of strategy 1
ant2 <- seq(0, 1, 1/nant) # Antigenicity of strategy 2
protection <- seq(0, 1, 1/nstrat) # "Strengths" of immune blockade
suppression <- seq(0, 1, 1/nstrat) # "Strengths" of immune suppression

param_grid <- meshgridn(list(ant1, protection, suppression))

#### Determine the maximum size of each strategy. Figure 1A ####
max_sizes_df <- data.frame("Antigenicity"=param_grid[[1]], "Blockade"=param_grid[[2]], "Suppression"=param_grid[[3]])
max_sizes_df$Max_size <- calc_max_size_vec(max_sizes_df$Antigenicity, max_sizes_df$Blockade, max_sizes_df$Suppression)
max_sizes_df$color <- rgb(0, 0, 0)
max_sizes_df$color[max_sizes_df$Max_size <= 0] <-pal_lancet()(8)[2]


max_sizes_df <- max_sizes_df %>% rename("S"= Suppression, "B"= Blockade)
max_sizes_p <- ggplot(max_sizes_df, aes(x=Antigenicity, y=Max_size, color=color)) +
  geom_segment(data = max_sizes_df, aes(x=Antigenicity, y=0, xend=Antigenicity, yend=Max_size)) +
  geom_point() +
  scale_color_identity() +
  ylab("Size") +
  ggtitle("Max Size") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text = element_text(face="bold", size=facet_font_size)) +
  coord_fixed(ratio=max(max_sizes_df$Antigenicity)/max(max_sizes_df$Max_size)) +
  guides(color=F)

max_sizes_p <- facet(max_sizes_p, facet.by = c("B", "S"),
                 short.panel.labs = FALSE)

ggsave("max_sizes.png", max_sizes_p, width = 10, height=10)

#### Get Lucky Vs Get Smart. Figure 1B ####
gl_vs_gs_sweep_grid <- meshgridn(list(ant1, ant2, protection, suppression))
gl_v_gs_df <- do.call(rbind,
                        lapply(seq_len(length(gl_vs_gs_sweep_grid[[1]])), function(i){
                          compete(p1_a = gl_vs_gs_sweep_grid[[2]][i],
                                  p1_p = gl_vs_gs_sweep_grid[[3]][i],
                                  p1_s = gl_vs_gs_sweep_grid[[4]][i],
                                  p1_name = "Get Smart",
                                  p2_a = gl_vs_gs_sweep_grid[[1]][i],
                                  p2_p = 0,
                                  p2_s = 0,
                                  p2_name = "Get Lucky")
                          }
                          )
                        )

gl_v_gs_df$result <- recode(gl_v_gs_df$result, "Get Lucky" = "GL", "Get Smart"="GS")


gl_v_gs_df <- gl_v_gs_df %>% rename("S"= "Get Smart_Sup", "B"= "Get Smart_Pro")
gl_v_gs_df$result[which(gl_v_gs_df$S==0 & gl_v_gs_df$B==0 & gl_v_gs_df$result %in% c("Co", "GS"))] <- "GL"
gl_v_gs_df$result <- factor(gl_v_gs_df$result, ordered = T, levels = c("GL", "GS", "Co", "N", "X"))

colnames(gl_v_gs_df) <- gsub(" ", "_", colnames(gl_v_gs_df))
gl_v_gs <- ggplot(gl_v_gs_df, aes(x=Get_Smart_Ant, y=Get_Lucky_Ant, fill=result)) +
  geom_tile() +
  xlab(paste(gsub("_", " ", unique(gl_v_gs_df$p1_name)), "Antigenicity")) +
  ylab(paste(gsub("_", " ", unique(gl_v_gs_df$p2_name)), "Antigenicity")) +
  ggtitle("Get Smart Vs Get Lucky") +
  scale_fill_lancet(name=NULL) +
  theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text = element_text(face="bold", size=facet_font_size)) +
  coord_equal()

gl_v_gs <- facet(gl_v_gs, facet.by = c( "B", "S"),
                     short.panel.labs = FALSE)

ggsave("gl_v_gs_comp_results.png", gl_v_gs, width = 12, height=12)


gl_v_gs_counts <- table(gl_v_gs_df$result)
gl_v_gs_counts <- gl_v_gs_counts/sum(gl_v_gs_counts)

gl_v_gs_counts <- ggplot(gl_v_gs_df, aes(x= result, fill=result)) +
  geom_histogram(stat="count") +
  ggtitle("Get Smart Vs Get Lucky") +
  scale_fill_lancet() +
  xlab(" ") +
  ylab("Count") +
  guides(fill=F)

#### Blockade Vs Suppression. Figure 1C ####
p_vs_s_sweep_grid <- meshgridn(list(ant1, ant2, suppression, protection))
p_v_s_df <- do.call(rbind,
                      lapply(seq_len(length(p_vs_s_sweep_grid[[1]])), function(i){
                        compete(p1_a = p_vs_s_sweep_grid[[2]][i],
                                p1_p = 0,
                                p1_s = p_vs_s_sweep_grid[[4]][i],
                                p1_name = "Suppression",
                                p2_a = p_vs_s_sweep_grid[[1]][i],
                                p2_p = p_vs_s_sweep_grid[[3]][i],
                                p2_s = 0,
                                p2_name = "Blockade")
                      }
                      )
)

p_v_s_df$result <- recode(p_v_s_df$result, "Suppression" = "S", "Blockade"="B")
p_v_s_df <- p_v_s_df %>% rename("S"= Suppression_Sup, "B"= Blockade_Pro)
p_v_s_df$result <- factor(p_v_s_df$result, ordered = T, levels = c("B", "S", "Co", "N", "X"))
p_v_s <- ggplot(p_v_s_df, aes(x= Suppression_Ant, y=Blockade_Ant, fill=result)) +
  geom_tile() +
  xlab(paste(gsub("_", " ", unique(p_v_s_df$p1_name)), "Antigenicity")) +
  ylab(paste(gsub("_", " ", unique(p_v_s_df$p2_name)), "Antigenicity")) +
  ggtitle("Suppression Vs Blockade") +
  scale_fill_lancet(name=NULL) +
  theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text = element_text(face="bold", size=facet_font_size)
        ) +
  coord_equal()

p_v_s <- facet(p_v_s, facet.by = c("B", "S"),
                 short.panel.labs = FALSE)


p_v_s_counts <- table(p_v_s_df$result)
p_v_s_counts <- p_v_s_counts/sum(p_v_s_counts)
signif(p_v_s_counts, 2)*100

p_v_s_counts <- ggplot(p_v_s_df, aes(x= result, fill=result)) +
  geom_histogram(stat="count") +
  ggtitle("Suppression Vs Blockade") +
  scale_fill_lancet() +
  xlab(" ") +
  ylab("Count") +
  guides(fill=F)

ggsave("p_v_s_comp_results.png", p_v_s, width = 12, height=12)
