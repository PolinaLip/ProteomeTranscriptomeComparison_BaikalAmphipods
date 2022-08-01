### Figure with up and down regulated proteins
library(cowplot)
library(ggplot2)
library(grid)
dir <- "~/labeglo2/proteome_transcr_comparision/slice_min_pvalue/"

species <- 'Gla'
data_up <- read.table(file.path(dir, 
           paste(species, '_AllupProteins_joinedWithTranscripts_sliceMinPadj.csv'))) # from signDEprot_with_RNAseqData.R
data_down <- read.table(file.path(dir, 
           paste(species, '_AlldownProteins_joinedWithTranscripts_sliceMinPadj.csv')))

data_up$method <- ifelse(data_up$method == 'MS/MS', 'Protein', 'Transcript')
data_down$method <- ifelse(data_down$method == 'MS/MS', 'Protein', 'Transcript')
data_up$method <- factor(data_up$method, 
                         levels = c('Transcript', 'Protein'))

data_up[data_up$condition == '24.6 °C',]$condition <- '24.7 °C'
data_up$condition <- factor(data_up$condition,
                            levels = c('6 °C', '24.7 °C'), 
                            labels = c('6 °C (control)', '24.7 °C (heat stress)'))
data_up$all_labels2 <- factor(data_up$all_labels,
                              levels=unique(
                                data_up$all_labels[order(data_up$pretty_varname)]))

data_down$method <- factor(data_down$method, 
                           levels = c('Transcript', 'Protein'))
data_down[data_down$condition == '24.6 °C',]$condition <- '24.7 °C'
data_down$condition <- factor(data_down$condition, 
                              levels = c('6 °C', '24.7 °C'),
                              labels = c('6 °C (control)', '24.7 °C (heat stress)'))
data_down$all_labels2 <- factor(data_down$all_labels,
                                levels=unique(
                                  data_down$all_labels[order(data_down$pretty_varname)]))

f <- function(x) {
  sapply(strsplit(x, '|', fixed=T), `[`, 2)
}
(gg_up <- ggplot(data_up, aes(x = method, y = values, color = condition)) +
    geom_point(position=position_jitterdodge(dodge.width=1),
               size = 0.7) +
    geom_boxplot(aes(fill = condition, 
                     linetype = sign2plot), 
                 outlier.alpha = 0, alpha = 0.4) +
    facet_wrap(~all_labels2, labeller = as_labeller(f), ncol = 5) +
    scale_linetype('adj. p-value:') +
    scale_color_manual('Condition:', values = c('#0571b0', '#ca0020')) +
    scale_fill_manual('Condition:', values = c('#0571b0', '#ca0020')) +
    theme_bw() +
    #ylab('Scaled absolute values') +
    ylab('Scaled absolute\nvalues') + #Gla
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12),
          strip.text = element_text(size = 8, margin = margin(2,2,2,2)))
)
(gg_down <- ggplot(data_down, aes(x = method, y = values, color = condition)) +
    geom_point(position=position_jitterdodge(dodge.width=1),
               size = 0.7) +
    geom_boxplot(aes(fill = condition, 
                     linetype = sign2plot), 
                 outlier.alpha = 0, alpha = 0.4) +
    facet_wrap(~all_labels2, labeller = as_labeller(f), ncol = 5) +
    scale_linetype('adj. p-value:') +
    scale_color_manual('Condition:', values = c('#0571b0', '#ca0020')) +
    scale_fill_manual('Condition:', values = c('#0571b0', '#ca0020')) +
    theme_bw() +
    #ylab('Scaled absolute values') +
    ylab('Scaled absolute\nvalues') + # Gla
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12),
          strip.text = element_text(size = 8, margin = margin(2,2,2,2)))
)

(together <- plot_grid(gg_up + theme(legend.position="none"),
                       #gg_down + theme(legend.position="none"), # Eve, Ecy
                       plot_grid(gg_down + theme(legend.position="none"), # Gla
                                 grid.rect(gp=gpar(col='white')), # Gla
                                 nrow=1, rel_widths = c(1, 1.23)), # Gla
                       labels = c("A", "B"),
                       label_y = 1.01,
                       hjust = -1,
                       nrow = 2
                       ,
                       rel_heights = c(0.9, 1) #Ecy, Gla
)
)
legend_b <- get_legend(gg_up + theme(legend.position="bottom", 
                                     legend.margin = margin(0,0,1,0)))

(together2 <- plot_grid(together, legend_b, ncol = 1, 
                        rel_heights = c(1, .05)))

ggsave(file.path(dir, paste(species, '_upAndDown_proteins_together_updated20220502.png')), 
       scale = 0.9, width = 9, height = 12) # Eve
ggsave(file.path(dir, paste(species, '_upAndDown_proteins_together_updated20220502.png')), 
       scale = 0.9, width = 9, height = 6.5) # Ecy
ggsave(file.path(dir, paste(species, '_upAndDown_proteins_together_updated20220502.png')), 
       scale = 0.9, width = 9, height = 4) # Gla

