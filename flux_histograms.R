#'
#'This R script plots the histograms of flux of all cell types per module
#'
#'
#'Input: flux file, barcodes and cell types, human module info file
#'
#'Output: histograms of flux of all cell types per module
#'
#'





rm(list = ls())
gc()

library(data.table)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(reshape)
library(ggridges)
library(ggplot2)


# load predicted flux---------------------------------------------------------
#COLON
path = "CLL_rendeiro_5_D0_flux.csv"
data_c <- read.csv(path)
#removes first column
data_c0 <- as.matrix(data_c[,-1])
#first column of data_c becomes row names of data_c0
rownames(data_c0) <- as.character(data_c[,1])
#transpose
data_c0 <- t(data_c0)
col_names_to_be_corrected <- colnames(data_c0)
col_names_corrected <- gsub('\\.', '-', col_names_to_be_corrected)
colnames(data_c0) <- col_names_corrected
ppp_all <-c()

#load factors with barcodes and cell types
load('CLL5_d0_ident_CYLD.RData') #CLL #5 day 0

#keep the elements of cell_id from column names of data_c0
yyy <- Idents_CLL5_D0[colnames(data_c0)]
for(ii in 1:nrow(data_c0)){
    xxx <- data_c0[ii,]
    final_df <- cbind(paste('X', 1:length(xxx), sep=''), xxx, yyy)
    final_df <- as.data.frame(final_df)
    final_df[,2] <- as.numeric(final_df[,2])
    colnames(final_df) <- c('var', 'flux', 'cellType')
    pp <- sd(final_df$flux)/abs(mean(final_df$flux))
    ppp_all <- c(ppp_all, pp)
}
#take those modules with npn zero flux (here set as 1e-10, re-adjust?)
tg_ids <- which(ppp_all > 1e-10)

# load mouse module info
load('mouse_module_info.RData')


# load human module info
human = read.csv('Human_M168_information.symbols.csv')
human$M_name = paste0(human$Compound_IN_name, '_',human$Compound_OUT_name) 
human_plot = human[, c(1, 2, 8, 3,4,5,6,7)]
human_plot = as.matrix(human_plot)
rownames(human_plot) = human_plot[,1]
human_plot = human_plot[,-1]
colnames(human_plot) = colnames(mouse_module_info)

#create list of modules number to be plotted e.g. 1 for M_1
plot_modules = c(61, 62, 63, 64, 65, 66, 67, 68)
plot_path_save = sub("_flux.csv.*", "_", path)  

for (ss in plot_modules){
    #insert module number of interest to get the row name position
    # jj is position/row of Module ss
    jj = which(tg_ids == ss) 
    #check if Module ss exists in tg_ids
    if(identical(jj, integer(0))){ 
        print(paste0("Module M_", ss," did not pass threshold."))
    }else{
        if(length(tg_ids) > 0){
            #jj = 20 # check module 2, runs on all modules
            xxx <- data_c0[paste0("M_",tg_ids[jj]), ]
            final_df <- data.frame(var = paste('X', 1:length(xxx), sep=''),
                                   flux = xxx,
                                   cellType = yyy)
            #human
            title <- human_plot[paste0("M_",tg_ids[jj]), 'M_name']
            
            aa <- ggplot(final_df, aes(x = flux, y = cellType, fill = cellType)) +
                geom_density_ridges() +
                theme_ridges() +
                theme(legend.position = 'none') +
                ggtitle(title) +
                theme(plot.title = element_text(hjust = 0.5))
            #plot(aa)
            ggsave(
                filename = paste0( plot_path_save ,"module_", ss,".pdf"),
                width = 13, height = 10, units = "in"
            )
        }
    }  
}