#### Figure with proteome/transcriptome correlations for all three species
library(cowplot)
library(gridExtra)
library(grid)
library(ggplot2)

p_format2 <- function(p) {
  if (p >= 0.05) {
    return(sprintf('p = %.3g', p))
  }
  rounded <- as.numeric(sprintf('%.1g', p))
  if (p >= rounded) {
    rounded <- rounded + 10 ^ floor(log10(rounded))
  }
  sprintf('p < %.1g', rounded)
}

dir_24 <- "~/labeglo2/proteome_transcr_comparision/slice_min_pvalue/"
dir_3 <- "~/labeglo2/proteome_transcr_comparision/slice_min_pvalue/"

x_min_limit <- -8.5 # for LFC with p-value < 0.05
x_min_limit <- -11.85 # for all observations

species_label_pos <- 2
stat_pos <- 1.3
label_size <- 5

### Eve, transcriptome from 24 hours ###
#eve_24h <- read.table(file.path(dir_24, 
#                                'Eve_24h_table_for_cor_plot_onlySign_sliceMinPValue.csv'))
eve_24h <- read.table(file.path(dir_24, 
                                'Eve_24h_table_for_cor_plot_All_sliceMinPValue.csv'))
cor_test_res_eve_24h <- cor.test(eve_24h$best_tlfc, 
                                 eve_24h$logFC)

max(eve_24h$best_tlfc)
min(eve_24h$best_tlfc)
(gg_eve_24 <- ggplot(eve_24h, aes(best_tlfc, logFC)) +
    geom_hline(yintercept = 0, color = '#387490', alpha = .7) +
    geom_vline(xintercept = 0, color = '#387490', alpha = .7) +
    geom_smooth(method = 'lm', color = 'grey55', fill = 'grey85') +
    geom_point(color='#009E73') +
    annotate(geom='text', x = x_min_limit, y = stat_pos, hjust = 0, 
             size = label_size,
             label = paste0('R = ', round(cor_test_res_eve_24h$estimate, 4), '\n',
                            p_format2(cor_test_res_eve_24h$p.value))) +
    annotate(geom='text', x = x_min_limit, y = species_label_pos, hjust = 0, 
             size = label_size,
             label = 'paste(italic(\"E.\"), \" \", italic(\"verrucosus\"))', 
             parse = T,
             color = '#009E73') +
    theme_light() +
    scale_x_continuous(limits = c(x_min_limit, 12.5)) +
    scale_y_continuous(limits = c(-1.2, 2.2)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_text(size = 14))
)

### Eve, transcriptome from 3 hours ###

#eve_3h <- read.table(file.path(dir_3, 
#                               'Eve_3h_table_for_cor_plot_onlySign_sliceMinPValue.csv'))
eve_3h <- read.table(file.path(dir_3, 'Eve_3h_table_for_cor_plot_All_sliceMinPValue.csv'))
cor_test_res_eve_3h <- cor.test(eve_3h$best_tlfc, 
                                eve_3h$logFC)

max(eve_3h$best_tlfc)
min(eve_3h$best_tlfc)
(gg_eve_3 <- ggplot(eve_3h, aes(best_tlfc, logFC)) +
    geom_hline(yintercept = 0, color = '#387490', alpha = .7) +
    geom_vline(xintercept = 0, color = '#387490', alpha = .7) +
    geom_smooth(method = 'lm', color = 'grey55', fill = 'grey85') +
    geom_point(color='#009E73') +
    annotate(geom='text', x = x_min_limit, y = stat_pos, hjust = 0, 
             size = label_size,
             label = paste0('R = ', round(cor_test_res_eve_3h$estimate, 4), '\n',
                            p_format2(cor_test_res_eve_3h$p.value))) + # 3h all proteins
    annotate(geom='text', x = x_min_limit, y = species_label_pos, hjust = 0, 
             size = label_size,
             label = 'paste(italic(\"E.\"), \" \", italic(\"verrucosus\"))', 
             parse = T,
             color = '#009E73') +
    theme_light() +
    scale_x_continuous(limits = c(x_min_limit, 12.5)) +
    scale_y_continuous(limits = c(-1.2, 2.2)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_text(size = 14))
)

### Ecy, transcriptome from 24 hours ###
#ecy_24h <- read.table(file.path(dir_24, 
#                                'Ecy_24h_table_for_cor_plot_onlySign_sliceMinPValue.csv'))
ecy_24h <- read.table(file.path(dir_24, 
                                'Ecy_24h_table_for_cor_plot_All_sliceMinPValue.csv'))
cor_test_res_ecy_24h <- cor.test(ecy_24h$best_tlfc, 
                                 ecy_24h$logFC)
max(ecy_24h$best_tlfc)
min(ecy_24h$best_tlfc)
(gg_ecy_24 <- ggplot(ecy_24h, aes(best_tlfc, logFC)) +
    geom_hline(yintercept = 0, color = '#387490', alpha = .7) +
    geom_vline(xintercept = 0, color = '#387490', alpha = .7) +
    geom_smooth(method = 'lm', color = 'grey55', fill = 'grey85') +
    #geom_point(color='#56B4E9') +
    geom_point(color='#148FD7') +
    annotate(geom='text', x = x_min_limit, y = stat_pos, hjust = 0, 
             size = label_size,
             label = paste0('R = ', round(cor_test_res_ecy_24h$estimate, 4), '\n',
                            p_format2(cor_test_res_ecy_24h$p.value))) + # 24h, all 
    annotate(geom='text', x = x_min_limit, y = species_label_pos, hjust = 0, 
             size = label_size,
             label = 'paste(italic(\"E.\"), \" \", italic(\"cyaneus\"))', 
             parse = T,
             color = '#148FD7') +
    theme_light() +
    scale_x_continuous(limits = c(x_min_limit, 12.5)) +
    scale_y_continuous(limits = c(-1.2, 2.2)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_text(size = 14))
)

### Ecy, transcriptome from 3 hours ###
#ecy_3h <- read.table(file.path(dir_3, 
#                               'Ecy_3h_table_for_cor_plot_onlySign_sliceMinPValue.csv'))
ecy_3h <- read.table(file.path(dir_3, 'Ecy_3h_table_for_cor_plot_All_sliceMinPValue.csv'))
cor_test_res_ecy_3h <- cor.test(ecy_3h$best_tlfc, 
                                ecy_3h$logFC)
max(ecy_3h$best_tlfc)
min(ecy_3h$best_tlfc)
(gg_ecy_3 <- ggplot(ecy_3h, aes(best_tlfc, logFC)) +
    geom_hline(yintercept = 0, color = '#387490', alpha = .7) +
    geom_vline(xintercept = 0, color = '#387490', alpha = .7) +
    geom_smooth(method = 'lm', color = 'grey55', fill = 'grey85') +
    geom_point(color='#148FD7') +
    annotate(geom='text', x = x_min_limit, y = stat_pos, hjust = 0, 
             size = label_size,
             label = paste0('R = ', round(cor_test_res_ecy_3h$estimate, 4), '\n',
                            p_format2(cor_test_res_ecy_3h$p.value))) + # 3h, all
    annotate(geom='text', x = x_min_limit, y = species_label_pos, hjust = 0, 
             size = label_size,
             label = 'paste(italic(\"E.\"), \" \", italic(\"cyaneus\"))', 
             parse = T,
             color = '#148FD7') +
    theme_light() +
    scale_x_continuous(limits = c(x_min_limit, 12.5)) +
    scale_y_continuous(limits = c(-1.2, 2.2)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_text(size = 14))
)
### Gla, transcriptome from 24 hours ###
#gla_24h <- read.table(file.path(dir_24, 
#                                'Gla_24h_table_for_cor_plot_onlySign_sliceMinPValue.csv'))
gla_24h <- read.table(file.path(dir_24, 'Gla_24h_table_for_cor_plot_All_sliceMinPValue.csv'))
cor_test_res_gla_24h <- cor.test(gla_24h$best_tlfc, 
                                 gla_24h$logFC)

max(gla_24h$best_tlfc)
min(gla_24h$best_tlfc)
(gg_gla_24 <- ggplot(gla_24h, aes(best_tlfc, logFC)) +
    geom_hline(yintercept = 0, color = '#387490', alpha = .7) +
    geom_vline(xintercept = 0, color = '#387490', alpha = .7) +
    geom_smooth(method = 'lm', color = 'grey55', fill = 'grey85') +
    geom_point(color='#E69F00') +
    annotate(geom='text', x = x_min_limit, y = stat_pos, hjust = 0, 
             size = label_size,
             label = paste0('R = ', round(cor_test_res_gla_24h$estimate, 4), '\n',
                            p_format2(cor_test_res_gla_24h$p.value))) +
    annotate(geom='text', x = x_min_limit, y = species_label_pos, hjust = 0, 
             size = label_size,
             label = 'paste(italic(\"G.\"), \" \", italic(\"lacustris\"))', 
             parse = T,
             color = '#E69F00') +
    theme_light() +
    scale_x_continuous(limits = c(x_min_limit, 12.5)) +
    scale_y_continuous(limits = c(-1.2, 2.2)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_text(size = 14))
)

### Gla, transcriptome from 3 hours ###
#gla_3h <- read.table(file.path(dir_3, 
#                               'Gla_3h_table_for_cor_plot_onlySign_sliceMinPValue.csv'))
gla_3h <- read.table(file.path(dir_3, 'Gla_3h_table_for_cor_plot_All_sliceMinPValue.csv'))
cor_test_res_gla_3h <- cor.test(gla_3h$best_tlfc, 
                                gla_3h$logFC)

max(gla_3h$best_tlfc)
min(gla_3h$best_tlfc)
(gg_gla_3 <- ggplot(gla_3h, aes(best_tlfc, logFC)) +
    geom_hline(yintercept = 0, color = '#387490', alpha = .7) +
    geom_vline(xintercept = 0, color = '#387490', alpha = .7) +
    geom_smooth(method = 'lm', color = 'grey55', fill = 'grey85') +
    geom_point(color='#E69F00') +
    annotate(geom='text', x = x_min_limit, y = stat_pos, hjust = 0, 
             size = label_size,
             label = paste0('R = ', round(cor_test_res_gla_3h$estimate, 4), '\n',
                            p_format2(cor_test_res_gla_3h$p.value))) +
    annotate(geom='text', x = x_min_limit, y = species_label_pos, hjust = 0,
             size = label_size,
             label = 'paste(italic(\"G.\"), \" \", italic(\"lacustris\"))', 
             parse = T,
             color = '#E69F00') +
    theme_light() +
    scale_x_continuous(limits = c(x_min_limit, 12.5)) +
    scale_y_continuous(limits = c(-1.2, 2.2)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_text(size = 14))
)

##### Put all plots together

add_axis <- function(g, x_name, y_name, ...) {
  x.grob <- textGrob(x_name, gp = gpar(...))
  y.grob <- textGrob(y_name, rot=90, gp = gpar(...))
  arrangeGrob(g, left = y.grob, bottom = x.grob)
}

g_top <- plot_grid(gg_ecy_24, gg_eve_24, gg_gla_24, nrow=1)
g_bot <- plot_grid(gg_ecy_3, gg_eve_3, gg_gla_3, nrow=1)

g_top <- plot_grid(gg_eve_24, gg_ecy_24, gg_gla_24, nrow=1) # different order: Eve-Ecy-Gla
g_bot <- plot_grid(gg_eve_3, gg_ecy_3, gg_gla_3, nrow=1) # different order: Eve-Ecy-Gla

plot_grid(add_axis(g_top, 'log2FC (24.6 °C/6 °C) transcriptome, 24 hours exposure', 
                   'log2FC (24.6 °C/6 °C)\nproteome, 24 hours exposure', fontsize = 17),
          grid.rect(gp=gpar(col=NA)),
          add_axis(g_bot, 'log2FC (24.6 °C/6 °C) transcriptome, 3 hours exposure', 
                   'log2FC (24.6 °C/6 °C)\nproteome, 24 hours exposure', fontsize = 17),
          ncol=1, labels = c('A', '', 'B'), label_size = 20,
          rel_heights = c(1.1, 0.05, 1.1), 
          label_y = 1.03)

(g_both <- plot_grid(
  grid.rect(gp=gpar(col=NA)),
  add_axis(g_top, 'log2FC (24.6 °C/6 °C) transcriptome, 24 hours exposure', 
           NULL, fontsize = 16),
  grid.rect(gp=gpar(col=NA)),
  add_axis(g_bot, 'log2FC (24.6 °C/6 °C) transcriptome, 3 hours exposure', 
           NULL, fontsize = 16),
  ncol=1, 
  labels = c('', 'A', '', 'B'), label_size = 20, label_x = 0, label_y = 1.11,
  rel_heights = c(0.1, 1.1, 0.05, 1.1)))
(g_both2 <- add_axis(g_both, NULL,
                     '     log2FC (24.6 °C/6 °C) proteome, 24 hours exposure',
                     fontsize = 16)); plot(g_both2)

ggsave(file.path(dir_24, 'all_species_cor_AllProteins_EveEcyGlaOrder.png'), g_both2,
       width = 15, height = 9, scale = 0.67)

#ggsave(file.path(dir_24, 'all_species_cor_pv05Proteins_EveEcyGlaOrder.png'), g_both2,
#       width = 15, height = 9, scale = 0.67)
