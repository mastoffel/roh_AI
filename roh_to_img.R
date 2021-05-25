library(tidyverse)
library(data.table)

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
roh_img_per_ind <- function(id) {
      
      roh_ind <- roh %>% filter(IID == id) 
      rohs <- unlist(Map(':', roh_ind$POS1_mb, roh_ind$POS2_mb))
      # make image
      pix[rohs] <- 255
      mat <- matrix(pix, nrow = side, byrow = T)
      
      surv <- juv_surv %>% 
                  dplyr::filter(id == {{ id }}) %>% 
                  .[["survival"]]
      
      # plot
      png(paste0("images/surv_", surv, "/", id, ".png")) 
      # 2. Create a plot
      par(mar=c(0,0,0,0))
      plot(mat, col = c("black", "white"), key = NULL, axis.col = NULL, axis.row=NULL,
           xlab = NA, ylab=NA, main=NA)
      # Close the pdf file
      dev.off() 
}

map(unique(juv_surv$id)[1:1000], roh_img_per_ind)
