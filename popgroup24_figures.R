library(dplyr)
library(magrittr)
library(pheatmap)
library(ggplot2)
library(ape)
library(ggtree)
library(reshape2)




##############################################
#          LTR distance histogram           #
############################################## 

# the fileswith dists here are made by script here:
# /scratch/botany/katie/te_ltrdist/get_ltr_dists.sh

impo_copia<-read.table("impo_copia_ltr_distance") %>% mutate(species="D. impolita", element="LTR Copia")
impo_gypsy<-read.table("impo_gypsy_ltr_distance") %>% mutate(species="D. impolita", element="LTR Gypsy")
sand_copia<-read.table("sand_copia_ltr_distance") %>% mutate(species="D. sandwicensis", element="LTR Copia")
sand_gypsy<-read.table("sand_gypsy_ltr_distance") %>% mutate(species="D. sandwicensis", element="LTR Gypsy")
yaho_copia<-read.table("yaho_copia_ltr_distance") %>% mutate(species="D. yahouensis", element="LTR Copia")
yaho_gypsy<-read.table("yaho_gypsy_ltr_distance") %>% mutate(species="D. yahouensis", element="LTR Gypsy")


all<-rbind(impo_copia, impo_gypsy, sand_copia, sand_gypsy, yaho_copia, yaho_gypsy)

all$species <- factor(all$species, levels=c("D. yahouensis", 
                                            "D. impolita", 
                                            "D. sandwicensis"))
png("popgroup24_ltr_distances.png", width = 1500, height = 1100)
ggplot(all, aes(x=V1, fill=species)) + 
  geom_histogram(position="identity", binwidth=0.001, alpha=0.5) +
  facet_wrap(~ element) +
  ylab("Count") +
  xlab("Percent Identity of 5' and 3' LTRs") +
  theme(text = element_text(size = 35), 
        legend.position = c(.18, .90),
        legend.title=element_blank(),
        legend.background=element_blank(),
        axis.text=element_text(size=30),
        axis.title=element_text(size=40),
        #axis.title.x=element_blank(),
        strip.text.x = element_text(size = 35),
        axis.text.x = element_text(angle = 45, size = 30, vjust = 1, hjust=1))
dev.off()





###############################################
#           draw orthofinder heatmap          #
############################################### 


# the orthogroups here were ge nerated on intact LTR sequences predicted by EDTA on yahoiuensis, impolita and sandwicensis
# orthofdinder was run with default params


orthologs<-read.table("/Users/katieemelianova/Desktop/Diospyros/presentations/popgroup24/Orthogroups.GeneCount.tsv", header=TRUE, col.names = c("D.impolita", "D.sandwicensis", "D.yahouensis", "total"))
orthologs %<>% dplyr::select("D.impolita", "D.sandwicensis", "D.yahouensis")

# get mean of TE orthogroup to scale with
ortho_sums<-orthologs %>% rowMeans()

colnames(orthologs) <- c("D. impolita", "D. sandwicensis", "D. yahouensis")
# draw heatmap scaling by< mean of orthogroup
png("popgroup24_orthogroup_heatmap.png", width = 1000, height = 900)
pheatmap((orthologs/ortho_sums), show_rownames=FALSE, treeheight_col=0, treeheight_row=0, fontsize_col=14, angle_col=45, border_color="black", cellwidth=200, cellheight=0.05)
dev.off()

library("RColorBrewer")
library(ComplexHeatmap)

ht = Heatmap(orthologs/ortho_sums, show_column_dend = FALSE, show_row_dend = FALSE, show_row_names = FALSE, show_heatmap_legend = FALSE, column_names_gp = grid::gpar(fontsize = 40), column_names_rot = 45)

library(circlize)
col_fun = colorRamp2(c(min(orthologs/ortho_sums), max(orthologs/ortho_sums)/2, max(orthologs/ortho_sums)), c("blue", "white", "red"))
lgd = Legend(col_fun = col_fun, title = "TE Copy Number", legend_height = unit(15, "cm"), grid_width = unit(3, "cm"), labels_gp = gpar(fontsize = 18), title_gp = gpar(fontsize = 30), labels = c("", "", ""), title_gap = unit(4, "mm"))
draw(lgd, x = unit(1, "npc"), y = unit(1, "npc"), just = c("right", "top"))

png("popgroup24_orthogroup_heatmap_test.png", width = 1300, height = 1200)
draw(ht, padding = unit(c(25, 30, 2, 110), "mm"))
draw(lgd, x = unit(0.99, "npc"), y = unit(0.99, "npc"), just = c("right", "top"))
dev.off()

######################################################
#       draw species specific TEs in orthogroups     #
######################################################

orthogroup_distances<-read.table("orthogroup_distances.tsv", col.names = c("species", "orthogroup", "TE", "start", "end", "distance"))

impolita_specific <- orthologs %>% filter(`D. impolita` > 1 & `D. yahouensis` < 1 & `D. sandwicensis` < 1)
yahouensis_specific<- orthologs %>% filter(`D. yahouensis` > 1 & `D. impolita` < 1 & `D. sandwicensis` < 1)
sandwicensis_specific <- orthologs %>% filter(`D. sandwicensis` > 1 & `D. impolita` < 1 & `D. yahouensis` < 1)
impolita_specific_cols<-impolita_specific %>% rownames()
yahouensis_specific_cols<-yahouensis_specific %>% rownames()
sandwicensis_specific_cols <- sandwicensis_specific %>% rownames()

impolita_specific_df<-orthogroup_distances %>% filter(orthogroup %in% impolita_specific_cols) %>% dplyr::select(species, distance) %>% mutate(cluster="D. impolita specific")
yahouensis_specific_df<-orthogroup_distances %>% filter(orthogroup %in% yahouensis_specific_cols) %>% dplyr::select(species, distance) %>% mutate(cluster="D. yahouensis specific")
sandwicensis_specific_df<-orthogroup_distances %>% filter(orthogroup %in% sandwicensis_specific_cols) %>% dplyr::select(species, distance) %>% mutate(cluster="D. sandwicensis specific")

species_specific_df<-rbind(impolita_specific_df, yahouensis_specific_df, sandwicensis_specific_df)

png("popgroup24_species_specific_dists.png", width = 800, height = 1000)
ggplot(species_specific_df, aes(x=distance, fill=cluster))+
  geom_histogram(alpha=0.5, position="dodge", bins = 90) +
  facet_wrap(~ cluster, ncol=1) +
  theme(legend.title=element_blank(),
        text = element_text(size = 30),
        #legend.position = c(.30, .9),
        legend.position = "none",
        strip.text.x = element_text(size = 30),) +
  ylab("Number of TEs") +
  xlab("Percent Identity of 5' and 3' LTRs") +
  scale_fill_manual(values=c("darkgreen", 
                             "dodgerblue", 
                             "firebrick1")) 
dev.off()


###################################################################
#         draw TEs in each species shared with sandwicensis       #
###################################################################


yahou_impo_shared<- orthologs %>% filter(`D. impolita` > 1 & `D. yahouensis` > 1 & `D. sandwicensis` < 1)
impo_sand_shared <- orthologs %>% filter(`D. impolita` > 1 &  `D. sandwicensis` > 1 & `D. yahouensis` < 1)
yahou_sand_shared <-orthologs %>% filter(`D. yahouensis` > 1 &  `D. sandwicensis` > 1 & `D. impolita` < 1)
yahou_impo_shared_cols<-yahou_impo_shared %>% rownames()
impo_sand_shared_cols<-impo_sand_shared %>% rownames()
yahou_sand_shared_cols<-yahou_sand_shared %>% rownames()
yahou_impo_shared_df<-orthogroup_distances %>% filter(orthogroup %in% yahou_impo_shared_cols) %>% dplyr::select(species, distance) %>% mutate(cluster="D. yahouensis & D. impolita shared")
impo_sand_shared_df<-orthogroup_distances %>% filter(orthogroup %in% impo_sand_shared_cols) %>% dplyr::select(species, distance) %>% mutate(cluster="D. sandwicensis & D. impolita shared")
yahou_sand_shared_df<-orthogroup_distances %>% filter(orthogroup %in% yahou_sand_shared_cols) %>% dplyr::select(species, distance) %>% mutate(cluster="D. yahouensis & D. sandwicensis shared")

species_shared_df<-rbind(impo_sand_shared_df, yahou_sand_shared_df)

png("popgroup24_sandwicensis_shared_dists.png", width = 800, height = 1000)
ggplot(species_shared_df, aes(x=distance, fill=cluster))+
  geom_histogram(alpha=0.5, position="dodge", bins = 90) +
  facet_wrap(~ cluster, ncol=1) +
  theme(legend.title=element_blank(),
        text = element_text(size = 40),
        #legend.position = c(.30, .9),
        legend.position = "none",
        strip.text.x = element_text(size = 30),) +
  ylab("Count") +
  xlab("Percent Identity of 5' and 3' LTRs")
dev.off()






yahou_shared<-yahou_impo_yellow_df %>% filter(species == "yahouensis") %>% mutate(cluster2="yahou_shared")
impo_shared<-yahou_impo_yellow_df %>% filter(species == "impolita") %>% mutate(cluster2="impo_shared")
impo_specific <- impo_red_df %>% mutate(cluster2="impo_specific")
yahou_specific <- yahou_red_df %>% mutate(cluster2="yahou_specific")

both_shared<-rbind(yahou_shared, impo_shared, impo_specific, yahou_specific)
both_shared[both_shared$species == "yahouensis",]$species = "D. yahouensis"
both_shared[both_shared$species == "impolita",]$species = "D. impolita"


png("popgroup24_impo_yahou_shared_unique_dists.png", width = 700, height = 900)
ggplot(both_shared, aes(x=distance, fill=cluster2))+
  geom_histogram(alpha=0.5, position="stack", bins = 90) +
  facet_wrap(~ species, ncol=1) +
  theme(legend.title=element_blank(),
        text = element_text(size = 30),
        #legend.position = c(.32, .9),
        legend.background=element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 35)) +
  ylab("Count")  +
  xlab("Percent Identity of 5' and 3' LTRs") +
  scale_fill_manual(labels= c("shared impolita and yahouensis", 
                              "specific to species", 
                              "shared", 
                              "specific"), 
                    values = c("hotpink", "darkgreen", "hotpink", "darkgreen"))
dev.off()




yahou_impo_yellow<- orthologs %>% filter(`D. impolita` > 1 & `D. yahouensis` > 1 & `D. sandwicensis` < 1)
imp_sand_yellow <- orthologs %>% filter(`D. impolita` > 1 &  `D. sandwicensis` > 1 & `D. yahouensis` < 1)
yahou_sand_yellow<-orthologs %>% filter(`D. yahouensis` > 1 &  `D. sandwicensis` > 1 & `D. impolita` < 1)
yahou_impo_yellow_cols<-yahou_impo_yellow %>% rownames()
imp_sand_yellow_cols<-imp_sand_yellow %>% rownames()
yahou_sand_yellow_cols<-yahou_sand_yellow %>% rownames()

yahou_impo_yellow_df<-orthogroup_distances %>% filter(orthogroup %in% yahou_impo_yellow_cols) %>% dplyr::select(species, distance) %>% mutate(cluster="D. yahouensis & D. impolita shared")
imp_sand_yellow_df<-orthogroup_distances %>% filter(orthogroup %in% imp_sand_yellow_cols) %>% dplyr::select(species, distance) %>% mutate(cluster="D. sandwicensis & D. impolita shared")
yahou_sand_yellow_df<-orthogroup_distances %>% filter(orthogroup %in% yahou_sand_yellow_cols) %>% dplyr::select(species, distance) %>% mutate(cluster="D. yahouensis & D. sandwicensis shared")


species_shared_df<-rbind(yahou_impo_yellow_df, imp_sand_yellow_df, yahou_sand_yellow_df)
species_shared_df<-rbind(imp_sand_yellow_df, yahou_sand_yellow_df)

png("popgroup24_species_shared_dists.png", width = 800, height = 1000)
ggplot(species_shared_df, aes(x=distance, fill=cluster))+
  geom_histogram(alpha=0.5, position="dodge", bins = 90) +
  facet_wrap(~ cluster, ncol=1) +
  theme(legend.title=element_blank(),
        text = element_text(size = 20),
        #legend.position = c(.30, .9),
        legend.position = "none",
        strip.text.x = element_text(size = 25),) +
  ylab("Count") +
  xlab("Percent Identity of 5' and 3' LTRs")
dev.off()


##################################################
#     barplot of TE numbers and genome size     #
################################################## 

# I here transcribe the numbers from the repeatmasker summary files stored in Diospyros/prsentations/EMBL*summary

san<-data.frame(species="D. sandwicensis",
                total_length=1007252419, 
                total_masked=726814292,
                Gypsy=402006824,
                Ty1Copia=94093519
                )

# subtract TE counts from categories which contain them
#san$total_length <-(san$total_length - (san$total_masked))
#san$total_masked <- (san$total_masked - (san$Ty1Copia + san$Gypsy))

imp<-data.frame(species="D. impolita",
                total_length=1834533962,
                total_masked=1388789270,
                Gypsy=833189331,
                Ty1Copia=224869633 
                )

yah<-data.frame(species="D. yahouensis",
                total_length=3015128523,
                total_masked=2273323726,
                Gypsy=1318885207,
                Ty1Copia=348328833
                )


# both both species together and reshape the df for plotting
# also convert units to Mb
both<-rbind(san, imp, yah) %>% melt()
both$value<-both$value/1000000

# orer species in order of size in ggplot
both$species <- factor(both$species, levels=c("D. yahouensis", 
                                            "D. impolita", 
                                            "D. sandwicensis"))


# plot barplot showing proportions of TEs in genome
png("popgroup24_TE_summary_barplot.png", width = 900, height = 900)
ggplot(both, aes(x=species, y=`value`, fill=variable)) + 
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(labels= c("Total Genome", 
                              "Total Masked",
                              "LTR Ty1 Copia", 
                              "LTR Gypsy"), 
                    values = c("royalblue1", "hotpink", "orange1", "palegreen3")) +
  theme(legend.title=element_blank(),
        axis.title.x=element_blank(),
        text = element_text(size = 50),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(size=40, angle=-25),
        legend.position = c(.80, .9)) +
  ylab("Total bases (Mbp)") 
dev.off()



##################################################
#     draw rooted highltighted species tree      #
################################################## 

species_tree<-ape::read.tree("/Users/katieemelianova/Desktop/Diospyros/diospyros_plots/lib1234_speciestree.nwk")
species_tree.rooted <- root(species_tree, which(species_tree$tip.label == "D.sandwicensis"))
tip<-c("D.olen", "D.fasciculosa", "D.macrocarpa", "D.ferrea")
species_tree.rooted<-drop.tip(species_tree.rooted, tip)
species_tree.rooted$species <- species_tree.rooted$tip.label


colours_tips <- case_when(species_tree.rooted$tip.label == "D.impolita" ~ "firebrick1",
                          species_tree.rooted$tip.label == "D.sandwicensis" ~"dodgerblue",
                          species_tree.rooted$tip.label == "D.yahouensis" ~"darkolivegreen3",
                          !(species_tree.rooted$tip.label %in% c("D.sandwicensis", "D.impolita", "D.yahouensis")) ~ "black")


dd <- data.frame(taxa=species_tree.rooted$tip.label, tipcols=colours_tips)
p<-ggtree(species_tree.rooted, size=3)
p <- p %<+% dd

png("popgroup24_species_highlight_tree.png", width = 1800, height = 1600)
p + geom_tiplab(size=22, aes(color=tipcols)) + 
  scale_colour_manual(values = c("black", "green4", "dodgerblue", "firebrick1")) +
  expand_limits(x = 0.07) + 
  #theme_tree2(plot.margin=margin(6, 120, 6, 6)) +
  theme(text=element_text(size=10), 
        legend.position = "none",
        #axis.line.x.bottom=element_line(size=3),
        #axis.text.x=element_text(size=40, vjust=-1),
        plot.margin = margin(, , 5, , "cm"))
  geom_treescale(x=0.06, y=3, fontsize=18, linesize=1)
dev.off()

# these refer to he plot.margin argument in theme above: which margins refer to which numbers
#up, #right, #down, #left



################################
#     plot soiltype tree       #
################################

# read in a modified version of the online spreadsheet; I removed some columns and renamed some of the soil types to match
# e.g. Seprpentines and Serpetine both become Serpentine
# update 22.08.23: we realised that a lot of the soil labels were incorrect and were taken from field notes. Correct mappings of soil to sample 
# are in the Paun et al. 2016 supplementary figure 5, though in a coloir coded tree form
# I took the designations from supp fig 5 and re-annotated the soil column based on this
# where there is no population data, soiltype is simply transposed. Where there is population data but all individuals in the supp figure 
# are the same, I transpose that soil type of all samples in the table. Where there are a mix of soiltypes per species
# and not all individuals are available in the tree for the table, not-included individuals in the table are left blank for soiltype
diospyros_localities<-read_tsv("/Users/katieemelianova/Desktop/Diospyros/diospyros_plots/Diospyros_Localities_corrected2016.tsv")

# Remove fields which dont have a location and dont have soiltype info
diospyros_localities %<>% filter(Latitude != 0 & !(is.na(Soil)))


# read in original tree
nwk<-ape::read.tree("/Users/katieemelianova/Desktop/Diospyros/diospyros_plots/RadiatingSpeciesDiospyros_ingrp-inds.nwk")

# split the tip name column by BT to get the population names
# I dont understand the format of this function and I dont like it at all
# but it works so whatevs
# also do some data cleaning
species<-sapply(strsplit(nwk$tip.label,"BT"), `[`, 1)
nwk$tip.label[nwk$tip.label =="impolitaBt103"] = "impolitaBT103"
sample_id<-sapply(strsplit(nwk$tip.label,"BT"), `[`, 2)
sample_id<-paste("BT", sample_id, sep = "")
sample_id[sample_id == "BTNA"] = nwk$tip.label[sample_id =="BTNA"]
sample_id<-sapply(strsplit(sample_id, "_"), `[`, 1)
sample_id<-sapply(strsplit(sample_id, "-"), `[`, 1)

# set the species in the newick object
nwk$species <- species
nwk$sample_id<-sample_id



# take the original data frame from the map part and take the rad and BT sample name columns only
mapping_loc<-diospyros_localities %>% 
  dplyr::select(`RAD localities`, `sequenced samples`, Soil) %>% 
  data.frame()

# this is a really dumb way to do it but whatever it works
# apply over each row, split the BTXXX column by ", ", and then add it to a data frame, where the other column is the RAD number
# this gives you a list of data frames, which are bound into one using bind_rows
bt_soil_mapping<-apply(mapping_loc, 1, function(x) strsplit(x[2], ", ") %>% data.frame(rad=x[3])) %>% bind_rows()

# remove the weird rownames the apply function gives it
rownames(bt_soil_mapping) <-NULL
colnames(bt_soil_mapping)<-c("sample_id", "soil")

# remove duplicate samples
bt_soil_mapping_noduplicates <- bt_soil_mapping %>% filter(sample_id != "BT147" & sample_id != "BT296")

# join tip names and soil types
bt_soil_mapping_joined<-left_join(data.frame(sample_id=nwk$sample_id), bt_soil_mapping_noduplicates) %>% dplyr::select(soil)

#set rownames to corresponding tip labels
rownames(bt_soil_mapping_joined) <- nwk$tip.label


# make a list where item name is species name and objects within are the tip labels belonging to that species
groupInfo<-split(nwk$tip.label, nwk$species)

# use groupOTU to group the tips by species and plot
nwk_grouped<-groupOTU(nwk, groupInfo, group_name = "species")

# take the same tree as above but make the tip labels the species name only (no population info)
nwk_grouped_speciesonly<-nwk_grouped
nwk_grouped_speciesonly$tip.label <- nwk_grouped_speciesonly$species


# first make one tree which is grouped and coloured by species
circ<-ggtree(nwk_grouped, aes(color=species), layout='circular') + geom_tiplab(size=5) 
circ_notiplab<-ggtree(nwk_grouped, aes(color=species), layout='circular', size=2.3) + geom_tiplab()
circ_notiplab<-ggtree(nwk_grouped, aes(color=species), layout='circular', size=2.3)


#offset=.012, width=.05
circ_soil<-gheatmap(circ_notiplab, 
                    bt_soil_mapping_joined, 
                    offset=0.0000005, width=.05, colnames=FALSE) + 
  scale_fill_viridis_d(option = "C", name = "Clade", na.value = "gray94") + 
  theme(legend.position="none")

# make the legend for the soiltypes
circ_nolegend<-ggtree(nwk)
soiltype_values<-bt_soil_mapping_joined %>% dplyr::select(soil)

soiltype_values$soil[rownames(soiltype_values) == "Ironcrust"] <- NA

# repeat this for soiltype
soiltype_legend<-gheatmap(circ_nolegend, soiltype_values, offset = .0001, width = 2,
                          colnames = FALSE,
                          colnames_offset_y = 1) +
  scale_fill_viridis_d(option = "C", na.value = "gray94", name = "soiltype", 
                       labels = c("Limestone", "Calcareous", "Schist", "Serpentine", "Ultramafic", "Volcanic", "Unknown")
                       #scale_fill_viridis_d(option = "C", na.value = "gray94", name = "Soiltype", labels = c("Ironcrust", "Kalcarious", "Serpentine", "Ultramafic", "Volcano-Sedimentary", "Unknown"))
                       #scale_fill_manual(values=mycolors, name="Soiltype"
  ) + theme(legend.text=element_text(size=50), 
            legend.title= element_blank())

# get the soiltype legend
soiltype_legend <- get_legend(soiltype_legend)

png("popgroup24_soiltype_tree.png", width = 1400, height = 1000)
ggarrange(circ_soil, soiltype_legend, widths = c(8,4))
dev.off()


