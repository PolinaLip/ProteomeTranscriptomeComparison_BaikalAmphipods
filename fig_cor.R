#### Figure with proteome/transcriptome correlation for all three species
library(cowplot)
library(gridExtra)
library(grid)

dir_24 <- "~/labeglo2/proteome_transcr_comparision"
dir_3 <- "~/labeglo2/proteome_transcr_comparision/3h/"

x_min_limit <- -8.5

# Eve, transcriptome from 24 hours
eve_24h <- read.table(file.path(dir_24, 'Eve_24h_table_for_cor_plot.csv'))
cor_test_res_eve_24h <- cor.test(eve_24h$best_tlfc, 
                                 eve_24h$logFC)

max(eve_24h$best_tlfc)
min(eve_24h$best_tlfc)
(gg_eve_24 <- ggplot(eve_24h, aes(best_tlfc, logFC)) +
  geom_hline(yintercept = 0, color = '#387490', alpha = .7) +
  geom_vline(xintercept = 0, color = '#387490', alpha = .7) +
  geom_smooth(method = 'lm', color = 'grey55', fill = 'grey85') +
  geom_point(color='#009E73') +
  annotate(geom='text', x = x_min_limit, y = 1.6, hjust = 0, size = 6,
           label = paste0('r2 = ', round(cor_test_res_eve_24h$estimate, 4), '\n',
                          'p-value ', p_format(cor_test_res_eve_24h$p.value,
                                               accuracy = 0.05))) +
  annotate(geom='text', x = x_min_limit, y = 2.2, hjust = 0, size = 6,
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

# Eve, transcriptome from 3 hours

eve_3h <- read.table(file.path(dir_3, 'Eve_3h_table_for_cor_plot.csv'))
cor_test_res_eve_3h <- cor.test(eve_3h$best_tlfc, 
                                 eve_3h$logFC)

max(eve_3h$best_tlfc)
min(eve_3h$best_tlfc)
(gg_eve_3 <- ggplot(eve_3h, aes(best_tlfc, logFC)) +
    geom_hline(yintercept = 0, color = '#387490', alpha = .7) +
    geom_vline(xintercept = 0, color = '#387490', alpha = .7) +
    geom_smooth(method = 'lm', color = 'grey55', fill = 'grey85') +
    geom_point(color='#009E73') +
    annotate(geom='text', x = x_min_limit, y = 1.6, hjust = 0, size = 6,
             label = paste0('r2 = ', round(cor_test_res_eve_3h$estimate, 4), '\n',
                            'p-value ', p_format(cor_test_res_eve_3h$p.value,
                                                 accuracy = 0.005))) +
    annotate(geom='text', x = x_min_limit, y = 2.2, hjust = 0, size = 6,
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

# Ecy, transcriptome from 24 hours
ecy_24h <- read.table(file.path(dir_24, 'Ecy_24h_table_for_cor_plot.csv'))
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
    annotate(geom='text', x = x_min_limit, y = 1.6, hjust = 0, size = 6,
             label = paste0('r2 = ', round(cor_test_res_ecy_24h$estimate, 4), '\n',
                            'p-value ', p_format(cor_test_res_ecy_24h$p.value,
                                                 accuracy = 0.00000001))) +
    annotate(geom='text', x = x_min_limit, y = 2.2, hjust = 0, size = 6,
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

# Ecy, transcriptome from 3 hours
ecy_3h <- read.table(file.path(dir_3, 'Ecy_3h_table_for_cor_plot.csv'))
cor_test_res_ecy_3h <- cor.test(ecy_3h$best_tlfc, 
                                ecy_3h$logFC)
max(ecy_3h$best_tlfc)
min(ecy_3h$best_tlfc)
(gg_ecy_3 <- ggplot(ecy_3h, aes(best_tlfc, logFC)) +
    geom_hline(yintercept = 0, color = '#387490', alpha = .7) +
    geom_vline(xintercept = 0, color = '#387490', alpha = .7) +
    geom_smooth(method = 'lm', color = 'grey55', fill = 'grey85') +
    geom_point(color='#148FD7') +
    annotate(geom='text', x = x_min_limit, y = 1.6, hjust = 0, size = 6,
             label = paste0('r2 = ', round(cor_test_res_ecy_3h$estimate, 4), '\n',
                            'p-value = ', p_format(cor_test_res_ecy_3h$p.value,
                                                 accuracy = 0.00000001))) +
    annotate(geom='text', x = x_min_limit, y = 2.2, hjust = 0, size = 6,
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
# Gla, transcriptome from 24 hours
gla_24h <- read.table(file.path(dir_24, 'Gla_24h_table_for_cor_plot.csv'))
cor_test_res_gla_24h <- cor.test(gla_24h$best_tlfc, 
                                gla_24h$logFC)

max(gla_24h$best_tlfc)
min(gla_24h$best_tlfc)
(gg_gla_24 <- ggplot(gla_24h, aes(best_tlfc, logFC)) +
    geom_hline(yintercept = 0, color = '#387490', alpha = .7) +
    geom_vline(xintercept = 0, color = '#387490', alpha = .7) +
    geom_smooth(method = 'lm', color = 'grey55', fill = 'grey85') +
    geom_point(color='#E69F00') +
    annotate(geom='text', x = x_min_limit, y = 1.6, hjust = 0, size = 6,
             label = paste0('r2 = ', round(cor_test_res_gla_24h$estimate, 4), '\n',
                            'p-value ', p_format(cor_test_res_gla_24h$p.value,
                                                 accuracy = 0.00001))) +
    annotate(geom='text', x = x_min_limit, y = 2.2, hjust = 0, size = 6,
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

# Gla, transcriptome from 3 hours
gla_3h <- read.table(file.path(dir_3, 'Gla_3h_table_for_cor_plot.csv'))
cor_test_res_gla_3h <- cor.test(gla_3h$best_tlfc, 
                                gla_3h$logFC)

max(gla_3h$best_tlfc)
min(gla_3h$best_tlfc)
(gg_gla_3 <- ggplot(gla_3h, aes(best_tlfc, logFC)) +
    geom_hline(yintercept = 0, color = '#387490', alpha = .7) +
    geom_vline(xintercept = 0, color = '#387490', alpha = .7) +
    geom_smooth(method = 'lm', color = 'grey55', fill = 'grey85') +
    geom_point(color='#E69F00') +
    annotate(geom='text', x = x_min_limit, y = 1.6, hjust = 0, size = 6,
             label = paste0('r2 = ', round(cor_test_res$estimate, 4), '\n',
                            'p-value ', p_format(cor_test_res$p.value,
                                                 accuracy = 0.05))) +
    annotate(geom='text', x = x_min_limit, y = 2.2, hjust = 0, size = 6,
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

#####

add_axis <- function(g, x_name, y_name, ...) {
  x.grob <- textGrob(x_name, gp = gpar(...))
  y.grob <- textGrob(y_name, rot=90, gp = gpar(...))
  arrangeGrob(g, left = y.grob, bottom = x.grob)
}

g_top <- plot_grid(gg_ecy_24, gg_eve_24, gg_gla_24, nrow=1)
g_bot <- plot_grid(gg_ecy_3, gg_eve_3, gg_gla_3, nrow=1)

plot_grid(add_axis(g_top, 'log2FC (24°C/6°C) transcriptome, 24 hours exposure', 
                   'log2FC (24°C/6°C) proteome,\n24 hours exposure', fontsize = 18),
          grid.rect(gp=gpar(col=NA)),
          add_axis(g_bot, 'log2FC (24°C/6°C) transcriptome, 3 hours exposure', 
                   'log2FC (24°C/6°C) proteome,\n24 hours exposure', fontsize = 18),
          ncol=1, labels = c('A', '', 'B'), label_size = 20,
          rel_heights = c(1, 0.05, 1), 
          label_y = 1.03)

ggsave(file.path(dir_24, 'all_species_cor_pv05.png'), scale = 2.8)


