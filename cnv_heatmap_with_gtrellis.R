#07/21/2022

#CNV heatmap for CNV output after analysis using 100Kb bin size
#The goals are:

# 1. convert the the log2 copy number gotten from CNVkit to absolute number
# 2. use heatmap to represent all the copy number change in the lumenal A samples
# 3. Then subtract or do copy ratio between primary and metastatic of each patient, to give you the CNV primary vs mets. 
# 4. Then plot these CNV for each sample using gtrellis too


library(dplyr)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
library(gtrellis)

#change working directory to where my files are 
setwd("/Users/olalekanusman/Desktop/PhD/Orgnoid_Project/troubleshooting_the_cnvkit_with_different_target_size/bin_size_100000")

#read in all the files and add a column that calculates the absolutes copy number.
pd5956a  <- read.csv("EGAN00001037943.cnr.txt", sep = "\t") %>% mutate(absolute_cn = 2 * (2**log2))
pd5956c  <- read.csv("EGAN00001270171.cnr.txt", sep = "\t") %>% mutate(absolute_cn = 2 * (2**log2))
pd9193a  <- read.csv("EGAN00001066968.cnr.txt", sep = "\t") %>% mutate(absolute_cn = 2 * (2**log2))
pd9193c  <- read.csv("EGAN00001270173.cnr.txt", sep = "\t") %>% mutate(absolute_cn = 2 * (2**log2))
pd9195a  <- read.csv("EGAN00001270176.cnr.txt", sep = "\t") %>% mutate(absolute_cn = 2 * (2**log2))
pd9195c  <- read.csv("EGAN00001270180.cnr.txt", sep = "\t") %>% mutate(absolute_cn = 2 * (2**log2))
pd11458a <- read.csv("EGAN00001270184.cnr.txt", sep = "\t") %>% mutate(absolute_cn = 2 * (2**log2))
pd11458c <- read.csv("EGAN00001270185.cnr.txt", sep = "\t") %>% mutate(absolute_cn = 2 * (2**log2))
pd11460a <- read.csv("EGAN00001270186.cnr.txt", sep = "\t") %>% mutate(absolute_cn = 2 * (2**log2))
pd11460c <- read.csv("EGAN00001270187.cnr.txt", sep = "\t") %>% mutate(absolute_cn = 2 * (2**log2))
pd11461a <- read.csv("EGAN00001270189.cnr.txt", sep = "\t") %>% mutate(absolute_cn = 2 * (2**log2))
pd11461c <- read.csv("EGAN00001270190.cnr.txt", sep = "\t") %>% mutate(absolute_cn = 2 * (2**log2))
pd11459a <- read.csv("EGAN00001270196.cnr.txt", sep = "\t") %>% mutate(absolute_cn = 2 * (2**log2))
pd11459c <- read.csv("EGAN00001270197.cnr.txt", sep = "\t") %>% mutate(absolute_cn = 2 * (2**log2))
pd13596a <- read.csv("EGAN00001270200.cnr.txt", sep = "\t") %>% mutate(absolute_cn = 2 * (2**log2))
pd13596c <- read.csv("EGAN00001270201.cnr.txt", sep = "\t") %>% mutate(absolute_cn = 2 * (2**log2))



col_fun = colorRamp2(c(0,1,2,3,4,5,6,8,10,12,14),c("#0571b0","#8cafd2","#f0f0f0","#fed5b9","#ffbb8b",
                                                   "#ffa061","#ff823a","#ff1e00","#c20704","#811006",
                                                   "#441104"))
#
col_fun = colorRamp2(c(0,1,2,3,4,5,6),c("#0571b0","#8cafd2","#f0f0f0","#ff823a","#ff1e00","#c20704","#811006"))


cm = ColorMapping(col_fun = col_fun)
lgd = color_mapping_legend(cm, plot = FALSE, title = "Chr copy")

gtrellis_layout(n_track = 17,
                track_axis = FALSE, 
                gap = unit(c(2,2),"mm"),
                track_height = unit.c(unit(6,"mm"),
                                      unit(1,"null"),unit(1,"null"),
                                      unit(1,"null"),unit(1,"null"),
                                      unit(1,"null"),unit(1,"null"),
                                      unit(1,"null"),unit(1,"null"),
                                      unit(1,"null"),unit(1,"null"),
                                      unit(1,"null"),unit(1,"null"),
                                      unit(1,"null"),unit(1,"null"),
                                      unit(1,"null"),unit(1,"null")),
                legend = lgd,
                add_ideogram_track = TRUE,
                ylab_rot = 0,
               
)

chr = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
all_chr = c(paste0("chr", 1:22),"chrX","chrY")
for(i in 1:length(chr)){
  add_track(category = all_chr[i], track = 1, panel_fun = function(gr) {
    grid.rect(gp = gpar(fill = "#EEEEEE"))
    grid.text(chr[i], gp = gpar(fontsize = 8))
  })
}

add_heatmap_track(pd5956a,pd5956a$absolute_cn,fill = col_fun)
add_heatmap_track(pd5956c,pd5956c$absolute_cn,fill = col_fun)
add_heatmap_track(pd9193a,pd9193a$absolute_cn,fill = col_fun)
add_heatmap_track(pd9193c,pd9193c$absolute_cn,fill = col_fun)
add_heatmap_track(pd9195a,pd9195a$absolute_cn,fill = col_fun)
add_heatmap_track(pd9195c,pd9195c$absolute_cn,fill = col_fun)
add_heatmap_track(pd11458a,pd11458a$absolute_cn,fill = col_fun)
add_heatmap_track(pd11458c,pd11458c$absolute_cn,fill = col_fun)
add_heatmap_track(pd11460a,pd11460a$absolute_cn,fill = col_fun)
add_heatmap_track(pd11460c,pd11460c$absolute_cn,fill = col_fun)
add_heatmap_track(pd11461a,pd11461a$absolute_cn,fill = col_fun)
add_heatmap_track(pd11461c,pd11461c$absolute_cn,fill = col_fun)
add_heatmap_track(pd11459a,pd11459a$absolute_cn,fill = col_fun)
add_heatmap_track(pd11459c,pd11459c$absolute_cn,fill = col_fun)
add_heatmap_track(pd13596a,pd13596a$absolute_cn,fill = col_fun)
add_heatmap_track(pd13596c,pd13596c$absolute_cn,fill = col_fun)


#Goal 3. calculating the fold change between the copy numnber of primary and met
pd5956 <- data.frame("chromosome"= pd5956a$chromosome,"start"= pd5956a$start, "end" = pd5956a$end,
                     "primary_cn"= pd5956a$absolute_cn,"met_cn"= pd5956c$absolute_cn) %>% 
  mutate(cn_fold_change = met_cn/primary_cn, log_cn_fold_change = log2(cn_fold_change))

pd9193 <- data.frame("chromosome"= pd9193a$chromosome,"start"= pd9193a$start, "end" = pd9193a$end,
                     "primary_cn"= pd9193a$absolute_cn,"met_cn"= pd9193c$absolute_cn) %>% 
  mutate(cn_fold_change = met_cn/primary_cn, log_cn_fold_change = log2(cn_fold_change))

pd9195 <- data.frame("chromosome"= pd9195a$chromosome,"start"= pd9195a$start, "end" = pd9195a$end,
                     "primary_cn"= pd9195a$absolute_cn,"met_cn"= pd9195c$absolute_cn) %>% 
  mutate(cn_fold_change = met_cn/primary_cn, log_cn_fold_change = log2(cn_fold_change))

pd11458 <- data.frame("chromosome"= pd11458a$chromosome,"start"= pd11458a$start, "end" = pd11458a$end,
                      "primary_cn"= pd11458a$absolute_cn,"met_cn"= pd11458c$absolute_cn) %>% 
  mutate(cn_fold_change = met_cn/primary_cn, log_cn_fold_change = log2(cn_fold_change))

pd11460 <- data.frame("chromosome"= pd11460a$chromosome,"start"= pd11460a$start, "end" = pd11460a$end,
                      "primary_cn"= pd11460a$absolute_cn,"met_cn"= pd11460c$absolute_cn) %>% 
  mutate(cn_fold_change = met_cn/primary_cn, log_cn_fold_change = log2(cn_fold_change))

pd11461 <- data.frame("chromosome"= pd11461a$chromosome,"start"= pd11461a$start, "end" = pd11461a$end,
                      "primary_cn"= pd11461a$absolute_cn,"met_cn"= pd11461c$absolute_cn) %>% 
  mutate(cn_fold_change = met_cn/primary_cn, log_cn_fold_change = log2(cn_fold_change))

pd11459 <- data.frame("chromosome"= pd11459a$chromosome,"start"= pd11459a$start, "end" = pd11459a$end,
                      "primary_cn"= pd11459a$absolute_cn,"met_cn"= pd11459c$absolute_cn) %>% 
  mutate(cn_fold_change = met_cn/primary_cn, log_cn_fold_change = log2(cn_fold_change))

pd13596 <- data.frame("chromosome"= pd13596a$chromosome,"start"= pd13596a$start, "end" = pd13596a$end,
                      "primary_cn"= pd13596a$absolute_cn,"met_cn"= pd13596c$absolute_cn) %>% 
  mutate(cn_fold_change = met_cn/primary_cn, log_cn_fold_change = log2(cn_fold_change))


#Goal 4. plotting heatmap for the fold change in chromosome between the primary and the met.

col_fun = colorRamp2(c(-1,0,1),c("#0571b0","#f0f0f0","#811006"))


cm = ColorMapping(col_fun = col_fun)
lgd = color_mapping_legend(cm, plot = FALSE, title = "log2(Fold Change)")

gtrellis_layout(n_track = 9,
                track_axis = FALSE, 
                gap = unit(2,"mm"),
                track_height = unit.c(unit(6,"mm"),
                                      unit(1,"null"),unit(1,"null"),
                                      unit(1,"null"),unit(1,"null"),
                                      unit(1,"null"),unit(1,"null"),
                                      unit(1,"null"),unit(1,"null")
                ),
                legend = lgd,
                axis_label_fontsize = 14,
                add_ideogram_track = TRUE
                
)

chr = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
all_chr = c(paste0("chr", 1:22),"chrX","chrY")
for(i in 1:length(chr)){
  add_track(category = all_chr[i], track = 1, panel_fun = function(gr) {
    grid.rect(gp = gpar(fill = "#EEEEEE"))
    grid.text(chr[i], gp = gpar(fontsize = 12))
  })
}

add_heatmap_track(pd5956,pd5956$log_cn_fold_change,fill = col_fun)
add_heatmap_track(pd9193,pd9193$log_cn_fold_change,fill = col_fun)
add_heatmap_track(pd9195,pd9195$log_cn_fold_change,fill = col_fun)
add_heatmap_track(pd11458,pd11458$log_cn_fold_change,fill = col_fun)
add_heatmap_track(pd11460,pd11460$log_cn_fold_change,fill = col_fun)
add_heatmap_track(pd11461,pd11461$log_cn_fold_change,fill = col_fun)
add_heatmap_track(pd11459,pd11459$log_cn_fold_change,fill = col_fun)
add_heatmap_track(pd13596,pd13596$log_cn_fold_change,fill = col_fun)




