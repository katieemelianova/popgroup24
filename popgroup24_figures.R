library(dplyr)
library(magrittr)

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

png("popgroup24_species_highlight_tree.png", width = 1500, height = 1500)
p + geom_tiplab(size=15, aes(color=tipcols)) + 
  scale_colour_manual(values = c("black", "darkolivegreen3", "dodgerblue", "firebrick1")) +
  expand_limits(x = 0.07) + 
  theme_tree2(plot.margin=margin(6, 120, 6, 6)) +
  theme(text=element_text(size=10), 
        legend.position = "none",
        axis.line.x.bottom=element_line(size=3),
        axis.text.x=element_text(size=40, vjust=-1),
        plot.margin = margin(, , 3, , "cm")) + 
  geom_treescale(x=0.06, y=3, fontsize=15, linesize=1)
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


