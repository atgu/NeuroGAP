##plotting up the 4 populations together so far: Moi, Mac, AAU, UCT

library(ggplot2)
library(gtable)
library(dplyr)
#install.packages('maptools')
library(maptools)
#install.packages('mapdata')#not installing. Download binary and install from there. timing out from CRAN
#instead trying rgdal as this may be the dependency now
library(mapdata) 
#library(rgdal)
#install.packages('raster')
library(raster)
#install.packages('ggmap')
library(ggmap)
#install.packages('cowplot')
library(cowplot)

#setwd('/Users/alicia/daly_lab/GINGER/agvp')
setwd("~/Dropbox/Projects/MGH-Broad/NeuroGap/pilotData/chipsCombined/")

pcs.AfricanPops <- read.table('/Users/elizabeth/Dropbox/Projects/MGH-Broad/NeuroGap/NeuroGAP/NeuroGAP_SharedWithZan/PCA/GSA_pilot/eigenvectors_tgp_agvp_neurogap_qc_ld_ancfilt', head=T)

pcs.NeuroGAP_AGVP_TGP <- read.table('/Users/elizabeth/Dropbox/Projects/MGH-Broad/NeuroGap/NeuroGAP/NeuroGAP_SharedWithZan/PCA/GSA_pilot/eigenvectors_tgp_neurogap_qc_ld', header=T)

colnames(pcs.AfricanPops) <- c("FID","IID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")


# Create Africa map -------------------------------------------------------
data(wrld_simpl)
world <- fortify(wrld_simpl)

# latitude/longitude bounds for africa
xlim = c(-20,50)
ylim = c(-35,35)

lims = SpatialPoints(coords = data_frame(x = xlim, y = ylim), proj4string = CRS("+proj=longlat +datum=WGS84"))%>%
  spTransform(CRS("+init=epsg:3857"))

# Read the latitude/longitdue data for populations in AGVP/TGP
pop_pos <- read.csv('african_pop_summary3_neurogap_combined.csv', header=T)#, stringsAsFactors = F)
#updated to change the Oromo in Ethopia symbol to an x to avoid confusion
pop_pos <- read.csv('pop_colors_neurogap2.csv', header=T)#, stringsAsFactors = F)
#pop_colors <- read.csv('pop_colors_neurogap2.csv', header=T)#, stringsAsFactors = F)

#pop_pos <- read.csv('pop_colors_neurogap_kemri.csv.csv', header=T)#, stringsAsFactors = F)
#pop_pos_plot <- pcs_labeled.AfricanPops
#pop_pos_plot <- subset(pcs_labeled.AfricanPops, Population != 'Ethiopia')

pop_pos_plot <- subset(pop_pos, Population != 'Ethiopia')
pop_pos_plot$Population <- factor(pop_pos_plot$Population, levels=as.character(pop_pos_plot$Population))
#pop_pos_plot <- subset(pop_pos_plot, Population %in% pca$FID)
color_vec <- as.character(pop_pos$Color)
names(color_vec) <- pop_pos$Population
shape_vec <- pop_pos$Shape
names(shape_vec) <- pop_pos$Population

# plot the map of Africa with data points labeled
#for the shapes defining language families, color_vec1 and shape_vec1
p1 <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group=group), fill='lightyellow', color='lightgrey') +
  geom_point(data = pop_pos_plot, aes(Longitude, Latitude, color=Population, fill=Population, shape=Population), size=2) +
  coord_fixed(xlim = c(-20,50), ylim = c(-35,35)) +
  labs(x='Longitude', y='Latitude') +
  theme_classic() +
  scale_fill_manual(name = "Population",
                    values = color_vec) +
  scale_color_manual(name = "Population",
                     values = color_vec) +
  scale_shape_manual(name = "Population",
                     values = shape_vec) +
  theme(panel.background = element_rect(fill = "lightblue"),
        legend.position='bottom')

p2 <- p1 + guides(fill=F, color=F, shape=F)

ggsave('AFR_AGVP_pos_allchips_updated.png', p1, width=8, height=5)

#don't bother with the language family encoding for the whole joint plot

# Plot PCA ----------------------------------------------------------------
#annotate the PCs with shape and color info we want, normal and with language shapes 
#pop_colors_neurogap <- pop_pos
#pop_colors_neurogap <- read.csv("pop_colors_neurogap2.csv", header=T)

pcs_labeled.Neuro <- merge(pcs.NeuroGAP_AGVP_TGP, pop_pos, by.x='FID', by.y='Population')
pcs_labeled.Neuro$Population <- factor(pcs_labeled.Neuro$FID)

#pop_colors_neurogap <- read.csv("african_pop_summary3_neurogap_uct.csv", header=T)

pcs_labeled.AfricanPops <- merge(pcs.AfricanPops, pop_colors_neurogap, by.x='FID', by.y='Population')
pcs_labeled.AfricanPops$Population <- factor(pcs_labeled.AfricanPops$FID)



# function for plotting PCA
#added a flag for filename to use for the PCs
pca_plot1 <- function(first_pc, second_pc, filename) {
  p_pca <- ggplot(filename, aes_string(x=first_pc, y=second_pc, color='FID', fill='FID', shape='FID')) +
    geom_point(alpha=1) +
    scale_color_manual(name = "Population",
                       values = color_vec, breaks="Population") +
    scale_fill_manual(name = "Population",
                      values = color_vec, breaks="Population") +
    scale_shape_manual(name = "Population",
                       values = shape_vec, breaks="Population") +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    theme_classic()
  return(p_pca)
}



#plot PCs 1 and 2 alone for the global subset
pc1_2_alone.Neuro <- ggplot(pcs_labeled.Neuro, aes_string(x='PC1', y='PC2', color='Population', fill='Population', shape='Population'), legendPosition="right") +
  geom_point(alpha=0.6) +
  scale_color_manual(name = "Population", breaks = pop_colors_neurogap$Population,values = color_vec) +
  scale_fill_manual(name = "Population", breaks = pop_colors_neurogap$Population, values = color_vec) +
  scale_shape_manual(name = "Population", breaks = pop_colors_neurogap$Population, values = shape_vec) +
  theme_classic()


#make sure looks right
#pca_plot1_2('PC1', 'PC2', pcs_labeled.Neuro) + guides(color=F, fill=F, shape=F)
pca_plot1('PC1', 'PC2', pcs_labeled.AfricanPops) + guides(color=F, fill=F, shape=F)
pca_plot1('PC1', 'PC2', pcs_labeled.Neuro) + guides(color=F, fill=F, shape=F)


# Make PCA biplots for the first 8 PCs
pc1_2 <- pca_plot1('PC1', 'PC2', pcs_labeled.Neuro) + guides(color=F, fill=F, shape=F)
pc3_4 <- pca_plot1('PC3', 'PC4', pcs_labeled.Neuro) + guides(color=F, fill=F, shape=F)
pc5_6 <- pca_plot1('PC5', 'PC6', pcs_labeled.Neuro) + guides(color=F, fill=F, shape=F)
pc7_8 <- pca_plot1('PC7', 'PC8', pcs_labeled.Neuro) + guides(color=F, fill=F, shape=F)

pc1_2.Africa <- pca_plot1('PC1', 'PC2', pcs_labeled.AfricanPops) + guides(color=F, fill=F, shape=F)
pc3_4.Africa <- pca_plot1('PC3', 'PC4', pcs_labeled.AfricanPops) + guides(color=F, fill=F, shape=F)
pc5_6.Africa <- pca_plot1('PC5', 'PC6', pcs_labeled.AfricanPops) + guides(color=F, fill=F, shape=F)
pc7_8.Africa <- pca_plot1('PC7', 'PC8', pcs_labeled.AfricanPops) + guides(color=F, fill=F, shape=F)


# Plot all of the PC together
pcs <- plot_grid(pc1_2, pc3_4, pc5_6, pc7_8, ncol=2, labels='AUTO', align='h')
pcs.Africa <- plot_grid(pc1_2.Africa, pc3_4.Africa, pc5_6.Africa, pc7_8.Africa, ncol=2, labels='AUTO', align='h')


# Draw a map with the PCA biplots
map_pcs <- ggdraw() +
  draw_plot(p1, x=0, y=0, width=0.5, height=1) +
  draw_plot(pcs, x=0.5, y=0, width=0.5, height=1)

map_pcs.Africa <- ggdraw() +
  draw_plot(p1, x=0, y=0, width=0.5, height=1) +
  draw_plot(pcs.Africa, x=0.5, y=0, width=0.5, height=1)


#and just the first 2 PCs blown up
map_pcs.Africa_1_2 <- ggdraw() +
  draw_plot(p1, x=0, y=0, width=0.5, height=1) +
  draw_plot(pc1_2.Africa, x=0.5, y=0, width=0.5, height=1)

map_pcs_1_2 <- ggdraw() +
  draw_plot(p1, x=0, y=0, width=0.5, height=1) +
  draw_plot(pc1_2, x=0.5, y=0, width=0.5, height=1)




#and for the global one without the map
pcs.1_2.Africa <- ggdraw() +
  draw_plot(pc1_2.Africa, x=0.5, y=0, width=0.5, height=1)

pcs.1_2 <- ggdraw() +
  draw_plot(pc1_2a, x=0.5, y=0, width=0.5, height=1)

# Save plots
save_plot('AFR_AGVP_NeuroGap_allplates_TGP_PCs.updated.pdf', pcs, base_height=7, base_width=8) #had to update the base to be wider to fit
save_plot('AFR1kg_AGVP_map_pca.allpilots.Africa.updated.pdf', map_pcs.Africa, base_height=7, base_width=14)

#save_plot('AFR1kg_AGVP_map_pca.UCT.Africa.pdf', map_pcs.Africa, base_height=7, base_width=14)


save_plot('AFR1kg_AGVP_map_pcs1_2.allpilots.Africa.pdf', map_pcs.Africa_1_2, base_height=7, base_width=14)
save_plot('AFR1kg_AGVP_map_pcs1_2.allpilots.TGP.updated.pdf', map_pcs_1_2, base_height=6, base_width=12)

pc1_2_alone.Neuro
