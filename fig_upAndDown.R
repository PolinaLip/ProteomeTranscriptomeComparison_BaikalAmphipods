### Figure with up and down regulated proteins
library(cowplot)
dir <- "~/labeglo2/proteome_transcr_comparision"

species <- 'Eve'
data_up <- read.table(file.path(dir, 
                                paste(species, '_AllupProteins_joinedWithTranscripts.csv')))
data_down <- read.table(file.path(dir, 
                                  paste(species, '_AlldownProteins_joinedWithTranscripts.csv')))

data_up$method <- factor(data_up$method, 
                               levels = c('RNAseq', 'MS/MS'))
data_up$condition <- factor(data_up$condition,
                                  levels = c('6 °C', '24.6 °C'))
data_up$all_labels2 <- factor(data_up$all_labels,
                              levels=unique(
                              data_up$all_labels[order(data_up$pretty_varname)]))

data_down$method <- factor(data_down$method, 
                         levels = c('RNAseq', 'MS/MS'))
data_down$condition <- factor(data_down$condition, 
                            levels = c('6 °C', '24.6 °C'))
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
  ylab('Scaled absolute values') +
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
  ylab('Scaled absolute values') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        strip.text = element_text(size = 8, margin = margin(2,2,2,2)))
)

(together <- plot_grid(gg_up + theme(legend.position="none"),
                   gg_down + theme(legend.position="none"),
                   #align = 'vh',
                   labels = c("A", "B"),
                   hjust = -1,
                   nrow = 2)
)
legend_b <- get_legend(gg_up + theme(legend.position="bottom"))

(together2 <- plot_grid(together, legend_b, ncol = 1, 
                        rel_heights = c(1, .05)))

ggsave(file.path(dir, paste(species, '_upAndDown_proteins_together.png')), 
       scale = 0.7, width = 12, height = 15)



