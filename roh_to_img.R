library(tidyverse)
library(data.table)
library(plot.matrix)
# fitness data
load("../sheep_roh/data/fitness_roh.RData") 
juv_surv <- fitness_data %>% 
                  filter(age == 0) %>% 
                  mutate(id = as.numeric(id)) %>% 
                  filter(!is.na(survival))

# genome_size
chrs <- read_delim("../sheep_ID/data/chromosome_info_oar31.txt", delim = "\t") %>% 
      rename(size_BP = Length,
             CHR = Part) %>% 
      mutate(size_KB = size_BP / 1000) %>% 
      filter(CHR %in% paste0("Chromosome ", 1:26)) %>% 
      mutate(CHR = as.numeric(str_replace(CHR, "Chromosome ", ""))) %>% 
      select(CHR, size_BP)

# cumulative bp positions
chrs_cum_length <- chrs %>% 
                     mutate(size_BP_cum = cumsum(lag(size_BP, default = 0)))

# roh
roh <- fread("../sheep_ID/output/ROH/roh.hom") %>% 
            left_join(chrs_cum_length) %>% 
            group_by(IID) %>% 
            mutate(POS1 = POS1 + size_BP_cum,
                   POS2 = POS2 + size_BP_cum,
                   POS1_mb = round(POS1/1e6, 0),
                   POS2_mb = round(POS2/1e6, 0)) 
      
# image size in pixels?           
genome <- sum(chrs$size_BP)/1e6
side <- round(sqrt(genome)+1, 0)
pix <- rep(0, side*side)

# save imgs in correct folder structure
roh_img_per_ind <- function(ind_id) {
      
      roh_ind <- roh %>% filter(IID == ind_id) 
      rohs <- unlist(Map(':', roh_ind$POS1_mb, roh_ind$POS2_mb))
      # make image
      pix[rohs] <- 255
      mat <- matrix(pix, nrow = side, byrow = T)
      
      surv <- juv_surv %>% 
                  dplyr::filter(id == ind_id) %>% 
                  .[["survival"]]
      
      # plot
      if (ind_id < 0) ind_id <- as.numeric(paste0("99", abs(ind_id)))
      
      png(paste0("images/surv_", surv, "/", ind_id, ".png")) 
      # 2. Create a plot
      par(mar=c(0,0,0,0))
      plot(mat, col = c("black", "white"), key = NULL, axis.col = NULL, axis.row=NULL,
           xlab = NA, ylab=NA, main=NA)
      # Close the pdf file
      dev.off() 
}

map(unique(juv_surv$id)[sample(1:nrow(juv_surv), 1000, replace=FALSE)], roh_img_per_ind)

l1 <- list.files("images/surv_0/")
l2 <- list.files("images/surv_1/")
