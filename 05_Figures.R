# Code used to generate figures in Manuscript
#Cleaned up Manuscript Figures
#
# (a) Preliminary Setup #####
library(ggplot2)
library(grid)
library(gridExtra)
library(MCMCvis)
library(mcmcplots)
library(coda)
library(dplyr)
library(tidyr)
library(broom)
library(ggspatial)
library(viridis)
library(cowplot)

path_to_data <- paste0(os_path, "Dropbox/Work Documents/occupancy_anacards_data/")
path_to_data2 <- "E:/"
memory.limit(1e13)

speciesCode <- "METTOX"
speciesName <- "Metopium toxiferum"

mod_cols <- c("#d01c8b", "#f1b6da", "#4dac26", "#b8e186")
mod_cols_full <- mod_cols[3:4]

load("E:/deduplicated_for_figures.RData")
deduplicated$citizenScience <- factor(deduplicated$citizenScience)
levels(deduplicated$citizenScience) <- c("Museum Specimen", "Citizen Science")
load(paste0(path_to_data, "Metopium toxiferum/surv_cov.RData"))
load(paste0(path_to_data, "florida.RData"))


width.cm <- 12.9
width.cm_onepanel<-8.4
width.cm_onepanel_small<-3.9
height.cm <- 6.45
pointsize <- 8

#(b) Create and save necessary dataframes for plotting #####

#* Collector covariates #####

collectors_df <- as.data.frame(sort(table(deduplicated$correctedCollector), decreasing=T))





# Fig 2: Collector Covariates #####




p1 <- ggplot(collectors_df, aes(x=Var1, y= Freq)) +
  geom_bar(stat="identity") +
  coord_flip() +
  ggtitle("(a) ") +
  theme(axis.text.y = element_blank()) +
  xlab("Collector") + ylab("Number of Records")+
  theme(text=element_text(size=10))



surv_cov$prev_detected <- factor(surv_cov$prev_detected)
levels(surv_cov$prev_detected) <- c("First Encounter", "Repeat Encounter")
levels(deduplicated$citizenScience) <- c("Herbarium Specimen", "Photograph")
temp <- subset(deduplicated, deduplicated$species == "Metopium toxiferum")
counts <- as.data.frame(table(temp$month))
decades <- as.data.frame(table(temp$decade))
citSci <- as.data.frame(table(temp$citizenScience))
detStat <- as.data.frame(table(surv_cov$prev_detected))



p2 <- ggplot(counts, aes(x=Var1, y= Freq)) +
  geom_bar(stat="identity") +
  ggtitle("(b) ") +
  xlab("Month") + ylab("Number of Records") +
  theme(text = element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90))

p3 <- ggplot(decades, aes(x=Var1, y= Freq)) +
  geom_bar(stat="identity") +
  ggtitle("(c) ") +
  xlab("Decade") + ylab("Number of Records") +
  theme(text = element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90))

p4 <- ggplot(citSci, aes(x=Var1, y= Freq)) +
  geom_bar(stat="identity") +
  ggtitle("(d) ") +
  xlab("Collection Type") + ylab("Number of Records") +
  theme(text = element_text(size=10))

p5 <- ggplot(detStat, aes(x=Var1, y= Freq)) +
  geom_bar(stat="identity") +
  ggtitle("(e) ") +
  xlab("Collection Type") + ylab("Number of Records") +
  theme(text = element_text(size=10))

abbott <- subset(deduplicated, deduplicated$correctedCollector =="J. R. Abbott")
#Abbott collected in 67 counties
lakela <- subset(deduplicated, deduplicated$correctedCollector =="O. K. Lakela")
#O. K. Lakela also collected in 34 counties
gillis <- subset(deduplicated, deduplicated$correctedCollector=="W. T. Gillis")

lak <- data.frame(county = unique(florida@data$NAME_2),
                  collected = rep(NA, 67))
abb <- data.frame(county=unique(florida@data$NAME_2),
                  collected = rep(NA, 67))
gil <- data.frame(county=unique(florida@data$NAME_2),
                  collected=rep(NA,67))
effort <- data.frame(county=unique(florida@data$NAME_2),
                     effort = as.vector(table(surv_cov$county)))

effort$log_effort <- log10(effort$effort)

for (i in 1:length(lak$county)) {
  lak$collected[i] <- lak$county[i] %in% unique(lakela$geo_county)
}

for (i in 1:length(abb$county)) {
  abb$collected[i] <- abb$county[i] %in% unique(abbott$geo_county)
}

for (i in 1:length(gil$county)) {
  gil$collected[i] <- gil$county[i] %in% unique(gillis$geo_county)
}



albersCRS <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'

floridaLakela <- sp::spTransform(florida, albersCRS)
floridaAbbott <- sp::spTransform(florida, albersCRS)
floridaGillis <- sp::spTransform(florida, albersCRS)
floridaEffort <- sp::spTransform(florida, albersCRS)

floridaLakela@data <- omnibus::insertCol(lak, into=florida@data, at='VARNAME_2', before=F)
floridaAbbott@data <- omnibus::insertCol(abb, into=florida@data, at='VARNAME_2', before=F)
floridaGillis@data <- omnibus::insertCol(gil, into=florida@data, at='VARNAME_2', before=F)
floridaEffort@data <- omnibus::insertCol(effort, into=florida@data, at='VARNAME_2', before=F)

id <- seq(1:67)
rownames(floridaLakela@data) <- id
id <- data.frame(id)
floridaLakela@data <- omnibus::insertCol(id, into=floridaLakela@data, at='collected', before=F)
Lakela_df <- tidy(floridaLakela, region = "id")
floridaLakela$id <- rownames(floridaLakela@data)
Lakela_df <- left_join(Lakela_df,
                       floridaLakela@data,
                       by = "id")

id <- seq(1:67)
rownames(floridaAbbott@data) <- id
id <- data.frame(id)
floridaAbbott@data <- omnibus::insertCol(id, into=floridaAbbott@data, at='collected', before=F)
Abbott_df <- tidy(floridaAbbott, region = "id")
floridaAbbott$id <- rownames(floridaAbbott@data)
Abbott_df <- left_join(Abbott_df,
                       floridaAbbott@data,
                       by = "id")

id <- seq(1:67)
rownames(floridaGillis@data) <- id
id <- data.frame(id)
floridaGillis@data <- omnibus::insertCol(id, into=floridaGillis@data, at='collected', before=F)
Gillis_df <- tidy(floridaGillis, region = "id")
floridaGillis$id <- rownames(floridaGillis@data)
Gillis_df <- left_join(Gillis_df,
                       floridaGillis@data,
                       by = "id")

id <- seq(1:67)
rownames(floridaEffort@data) <- id
id <- data.frame(id)
floridaEffort@data <- omnibus::insertCol(id, into=floridaEffort@data, at='log_effort', before=F)
Effort_df <- tidy(floridaEffort, region = "id")
floridaEffort$id <- rownames(floridaEffort@data)
Effort_df <- left_join(Effort_df,
                       floridaEffort@data,
                       by = "id")

p6 <- ggplot() +
  coord_sf(crs=albersCRS) +
  ggtitle("(f) ") +
  geom_polygon(data=Lakela_df, aes(x = long, y=lat, group=group, fill=collected)) +
  xlab("") + ylab("") +
  scale_fill_manual (values=c("white", "orange")) +
  geom_path(data = Lakela_df, aes(x = long, y = lat, group = group)) +
  annotation_scale(location = "bl", width_hint = 0.45, text_cex=.8) +
  annotation_north_arrow(location = "bl", which_north = "true",
                         pad_x = unit(0, "in"), pad_y = unit(0.4, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position="none") +
  theme(text = element_text(size=10))


p7 <- ggplot() +
  coord_sf(crs=albersCRS) +
  ggtitle("(g) ") +
  geom_polygon(data=Abbott_df, aes(x = long, y=lat, group=group, fill=collected)) +
  xlab("") + ylab("") +
  scale_fill_manual (values=c("white", "orange")) +
  geom_path(data = Lakela_df, aes(x = long, y = lat, group = group)) +
  annotation_scale(location = "bl", width_hint = 0.45, text_cex=.8) +
  annotation_north_arrow(location = "bl", which_north = "true",
                         pad_x = unit(0, "in"), pad_y = unit(0.4, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position="none") +
  theme(text = element_text(size=10))

p8_temp <- ggplot() +
  coord_sf(crs=albersCRS) +
  ggtitle("(h) ") +
  geom_polygon(data=Effort_df, aes(x = long, y=lat, group=group, fill=log_effort)) +
  xlab("") + ylab("") +
  labs(fill = "log (Number \n of Records)") +
  scale_fill_viridis() +
  geom_path(data = Effort_df, aes(x = long, y = lat, group = group)) +
  annotation_scale(location = "bl", width_hint = 0.45, text_cex=.8) +
  annotation_north_arrow(location = "bl", which_north = "true",
                         pad_x = unit(0, "in"), pad_y = unit(0.4, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(legend.title.align=0.5) +
  theme(text = element_text(size=10))

p8 <- p8_temp + theme(legend.position="none")

#Create a centered legend

g <- ggplotGrob(p8_temp + theme(legend.position="right"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

# extract guides table
guides_index <- which(sapply(legend$grobs, function(x) x$name) == "layout")
hjust <- 0.5
# there can be multiple guides within one legend box
for (gi in guides_index) {
  guides <- legend$grobs[[gi]]

  # add extra column for spacing
  # guides$width[5] is the extra spacing from the end of the legend text
  # to the end of the legend title. If we instead distribute it by `hjust:(1-hjust)` on
  # both sides, we get an aligned legend
  spacing <- guides$width[5]
  guides <- gtable::gtable_add_cols(guides, hjust*spacing, 1)
  guides$widths[6] <- (1-hjust)*spacing
  title_index <- guides$layout$name == "title"
  guides$layout$l[title_index] <- 2

  # reconstruct guides and write back
  legend$grobs[[gi]] <- guides
}





lay= rbind(c(1,1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3),
           c(1,1, 1,  4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5),
           c(1,1, 1, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9))


png("./collector_behavior2.png", width=8,
    height=10, units="in", pointsize=pointsize, res=300
)

grid.arrange(p1,p2,p3,p4,p5, p6, p7, p8,legend, layout_matrix=lay)
dev.off()

#Fig 3: Cross validation results ######


#MANIND
load(paste0(path_to_data, "MANIND/crossvalidation/MANIND_crossvalidation_results_summary.RData"))
x <- 1:4
MANIND_ll <- results_summary
MANIND_ll <- cbind(MANIND_ll, x)
mod_cols2 <- mod_cols[c(3,4,1,2)]

p_MANIND<-ggplot(data=MANIND_ll, aes(x=x, y=log_lik_scale, fill=method))+
  ylim(c(0,7)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=.5)) +
  geom_bar(stat="identity")+
  ylab("Scaled Log Score") + xlab("") +
  theme(text = element_text(size=12)) +
  scale_fill_manual(values=mod_cols2, name="Method",
                    labels=c("Full Model - \n All Data", "Full Model -  \n Filtered",
                             "Simple Model - \n All Data", "Simple Model -  \n Filtered")) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())

g <- ggplotGrob(p_MANIND + theme(legend.position="right"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]



#METTOX
load(paste0(path_to_data, "METTOX/crossvalidation/METTOX_crossvalidation_results_summary.RData"))
METTOX_ll <- results_summary

p_METTOX<-ggplot(data=METTOX_ll, aes(x=x, y=log_lik_scale, fill=method))+
  ylim(c(0,7)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=.5)) +
  geom_bar(stat="identity")+
  ylab("Scaled Log Score") + xlab("") +
  theme(text = element_text(size=20)) +
  scale_fill_manual(values=mod_cols2, name="Method",
                    labels=c("Full Model -  All Data", "Full Model -   Filtered",
                             "Simple Model -  All Data", "Simple Model -  Filtered")) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())

#RHUCOP
load(paste0(path_to_data, "RHUCOP/crossvalidation/RHUCOP_crossvalidation_results_summary.RData"))
RHUCOP_ll <- results_summary

p_RHUCOP<-ggplot(data=RHUCOP_ll, aes(x=x, y=log_lik_scale, fill=method))+
  ylim(c(0,7)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=.5)) +
  geom_bar(stat="identity")+
  ylab("Scaled Log Score") + xlab("") +
  theme(text = element_text(size=20)) +
  scale_fill_manual(values=mod_cols2, name="Method",
                    labels=c("Full Model -  All Data", "Full Model -   Filtered",
                             "Simple Model -  All Data", "Simple Model -  Filtered")) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
#SCHTER
load(paste0(path_to_data, "SCHTER/crossvalidation/SCHTER_crossvalidation_results_summary.RData"))
SCHTER_ll <- results_summary

p_SCHTER<-ggplot(data=SCHTER_ll, aes(x=x, y=log_lik_scale, fill=method))+
  ylim(c(0,7)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=.5)) +
  geom_bar(stat="identity")+
  ylab("Scaled Log Score") + xlab("") +
  theme(text = element_text(size=20)) +
  scale_fill_manual(values=mod_cols2, name="Method",
                    labels=c("Full Model -  All Data", "Full Model -   Filtered",
                             "Simple Model -  All Data", "Simple Model -  Filtered")) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
#TOXPUB
load(paste0(path_to_data, "TOXPUB/crossvalidation/TOXPUB_crossvalidation_results_summary.RData"))
TOXPUB_ll <- results_summary

p_TOXPUB<-ggplot(data=TOXPUB_ll, aes(x=x, y=log_lik_scale, fill=method))+
  ylim(c(0,7)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=.5)) +
  geom_bar(stat="identity")+
  ylab("Scaled Log Score") + xlab("") +
  theme(text = element_text(size=20)) +
  scale_fill_manual(values=mod_cols2, name="Method",
                    labels=c("Full Model -  All Data", "Full Model -   Filtered",
                             "Simple Model -  All Data", "Simple Model -  Filtered")) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
#TOXRAD
load(paste0(path_to_data, "TOXRAD/crossvalidation/TOXRAD_crossvalidation_results_summary.RData"))
TOXRAD_ll <- results_summary

p_TOXRAD<-ggplot(data=TOXRAD_ll, aes(x=x, y=log_lik_scale, fill=method))+
  ylim(c(0,7)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=.5)) +
  geom_bar(stat="identity")+
  ylab("Scaled Log Score") + xlab("") +
  theme(text = element_text(size=20)) +
  scale_fill_manual(values=mod_cols2, name="Method",
                    labels=c("Full Model -  All Data", "Full Model -   Filtered",
                             "Simple Model -  All Data", "Simple Model -  Filtered")) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
#TOXVER
load(paste0(path_to_data, "TOXVER/crossvalidation/TOXVER_crossvalidation_results_summary.RData"))
TOXVER_ll <- results_summary

p_TOXVER<-ggplot(data=TOXVER_ll, aes(x=x, y=log_lik_scale, fill=method))+
  ylim(c(0,7)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=.5)) +
  geom_bar(stat="identity")+
  ylab("Scaled Log Score") + xlab("") +
  theme(text = element_text(size=20)) +
  scale_fill_manual(values=mod_cols2, name="Method",
                    labels=c("Full Model -  All Data", "Full Model -   Filtered",
                             "Simple Model -  All Data", "Simple Model -  Filtered")) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())

min(MANIND_ll$log_lik_scale, METTOX_ll$log_lik_scale,
    RHUCOP_ll$log_lik_scale, SCHTER_ll$log_lik_scale,
    TOXPUB_ll$log_lik_scale, TOXRAD_ll$log_lik_scale,
    TOXVER_ll$log_lik_scale) #2.57

max(MANIND_ll$log_lik_scale, METTOX_ll$log_lik_scale,
    RHUCOP_ll$log_lik_scale, SCHTER_ll$log_lik_scale,
    TOXPUB_ll$log_lik_scale, TOXRAD_ll$log_lik_scale,
    TOXVER_ll$log_lik_scale) #6.82


p_MANIND2 <- p_MANIND + theme(legend.position = "none", text = element_text(size=12))
p_METTOX2 <- p_METTOX + theme(legend.position = "none", text = element_text(size=12))
p_RHUCOP2 <- p_RHUCOP + theme(legend.position = "none", text = element_text(size=12))
p_SCHTER2 <- p_SCHTER + theme(legend.position = "none", text = element_text(size=12))
p_TOXPUB2 <- p_TOXPUB + theme(legend.position = "none", text = element_text(size=12))
p_TOXRAD2 <- p_TOXRAD + theme(legend.position = "none", text = element_text(size=12))
p_TOXVER2 <- p_TOXVER + theme(legend.position = "none", text = element_text(size=12))

# png("./log_lik_barplot.png", width=7,
#     height=4, units="in", pointsize=pointsize, res=300
# )
# grid.arrange(p_MANIND,p_METTOX, p_RHUCOP, p_SCHTER, p_TOXPUB, p_TOXRAD, p_TOXVER, ncol=7,
#              widths=c(1, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75))
# dev.off()

grid.arrange(p_MANIND2, p_METTOX2, p_RHUCOP2, p_SCHTER2, p_TOXPUB2,
             p_TOXRAD2, p_TOXVER2, legend, ncol=4)
dev.off()


p_MANIND3 <- p_MANIND2 + theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) + ggtitle("(a)")
p_METTOX3 <- p_METTOX2 + theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.y=element_blank(),
                               axis.ticks.y=element_blank()) + ggtitle("(b)")
p_RHUCOP3 <- p_RHUCOP2 + theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.y=element_blank(),
                               axis.ticks.y=element_blank()) + ggtitle("(c)")
p_SCHTER3 <- p_SCHTER2 + theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.y=element_blank(),
                               axis.ticks.y=element_blank()) + ggtitle("(d)")
p_TOXPUB3 <- p_TOXPUB2 + theme(axis.title.x=element_blank(), axis.ticks.x = element_blank()) + ggtitle("(e)")
p_TOXRAD3 <- p_TOXRAD2 + theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.y=element_blank(),
                               axis.ticks.y=element_blank()) + ggtitle("(f)")
p_TOXVER3 <- p_TOXVER2 + theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.y=element_blank(),
                               axis.ticks.y=element_blank()) + ggtitle("(g)")


g_MANIND <- ggplotGrob(p_MANIND3)
g_METTOX <- ggplotGrob(p_METTOX3)
g_RHUCOP <- ggplotGrob(p_RHUCOP3)
g_SCHTER <- ggplotGrob(p_SCHTER3)
g_TOXPUB <- ggplotGrob(p_TOXPUB3)
g_TOXRAD <- ggplotGrob(p_TOXRAD3)
g_TOXVER <- ggplotGrob(p_TOXVER3)



# Set the widths
g_METTOX$widths <- g_MANIND$widths
g_RHUCOP$widths <- g_MANIND$widths
g_SCHTER$widths <- g_MANIND$widths
g_TOXPUB$widths <- g_MANIND$widths
g_TOXRAD$widths <- g_MANIND$widths
g_TOXVER$widths <- g_MANIND$widths

png("./log_lik_barplot.png", width=7, height=4, units="in", pointsize=pointsize, res=300)
grid.arrange(g_MANIND, g_METTOX, g_RHUCOP, g_SCHTER, g_TOXPUB, g_TOXRAD, g_TOXVER, legend, ncol=4)
dev.off()

#Fig 4: Detection Caterplots #####

speciesList <- c("MANIND", "METTOX", "RHUCOP", "SCHTER", "TOXPUB", "TOXRAD", "TOXVER")

df_zeta <- data.frame(species = rep(c("MANIND", "METTOX", "RHUCOP", "SCHTER", "TOXPUB", "TOXRAD", "TOXVER"), each= 2),
                      model = rep(c("full", "filter"), 7),
                      x = 1:14,
                      mean = rep(NA,14),
                      L1 = rep(NA, 14),
                      L2 = rep(NA, 14),
                      U2 = rep(NA,14),
                      U1 = rep(NA,14))

df_eta <- data.frame(species = rep(c("MANIND", "METTOX", "RHUCOP", "SCHTER", "TOXPUB", "TOXRAD", "TOXVER"), each= 2),
                     model = rep(c("full", "filter"), 7),
                     x = 1:14,
                     mean = rep(NA,14),
                     L1 = rep(NA, 14),
                     L2 = rep(NA, 14),
                     U2 = rep(NA,14),
                     U1 = rep(NA,14))


pos <- c(1, 3, 5, 7, 9, 11, 13)


for (i in 2:7) {
  load(paste0(path_to_data, speciesList[i], "/samples_1.RData"))
  chain1 <- samples1$samples
  rm(samples1)
  load(paste0(path_to_data, speciesList[i], "/samples_2.RData"))
  chain2 <- samples2$samples
  rm(samples2)
  mcmc <- mcmc.list(chain1, chain2)


  thing_all <- MCMCsummary(mcmc, params=c("zeta"),
                           probs=c(0.025, 0.16, 0.84, 0.975))[1, c(1, 3:6)]
  df_zeta[pos[i], 4:8] <- thing_all
  rm(mcmc, thing_all)
  load(paste0(path_to_data, speciesList[i], "/fullmod_filtered_mcmcList.RData"))
  thing_filter <- MCMCsummary(mcmc, params=c("zeta"),
                              probs=c(0.025, 0.16, 0.84, 0.975))[1, c(1, 3:6)]
  df_zeta[(pos[i]+1), 4:8] <- thing_filter
  rm(mcmc, thing_filter)
}

#Go back through for eta
for (i in 1:7) {
  load(paste0(path_to_data, speciesList[i], "/samples_1.RData"))
  chain1 <- samples1$samples
  rm(samples1)
  load(paste0(path_to_data, speciesList[i], "/samples_2.RData"))
  chain2 <- samples2$samples
  rm(samples2)
  mcmc <- mcmc.list(chain1, chain2)


  thing_all <- MCMCsummary(mcmc, params=c("eta"),
                           probs=c(0.025, 0.16, 0.84, 0.975))[1, c(1, 3:6)]
  df_eta[pos[i], 4:8] <- thing_all
  rm(mcmc, thing_all)
  load(paste0(path_to_data, speciesList[i], "/fullmod_filtered_mcmcList.RData"))
  thing_filter <- MCMCsummary(mcmc, params=c("eta"),
                              probs=c(0.025, 0.16, 0.84, 0.975))[1, c(1, 3:6)]
  df_eta[(pos[i]+1), 4:8] <- thing_filter
  rm(mcmc, thing_filter)
  cat(i)
}
min(df$L1) #1.26
max(df$U1) #4.97
min(df_eta$L1) #-4.99
max(df_eta$U1) #2.92

df_eta$inter <- NA
for(i in 1:length(df_eta$inter)){
  df_eta$inter[i] <- ( df_eta$L1[i]<0) && (df_eta$U1[i]>0)
}
df_eta$col <- c(mod_cols[3], mod_cols[4])
for(i in 1:length(df_eta$col)) {
  if(df_eta$inter[i]==1){df_eta$col[i] <- "grey"}
}

#zeta
df_zeta$inter <- NA
for(i in 1:length(df_zeta$inter)){
  df_zeta$inter[i] <- ( df_zeta$L1[i]<0) && (df_zeta$U1[i]>0)
}
df_zeta$col <- c(mod_cols[3], mod_cols[4])
for(i in 1:length(df_zeta$col)) {
  if(df_zeta$inter[i]==1){df_zeta$col[i] <- "grey"}
}


save(df_zeta, file=paste0(path_to_data, "crossval_zeta.RData"))
save(df_eta, file=paste0(path_to_data, "crossval_eta.RData"))

p1 <- ggplot(df_zeta, aes(x=x, y=mean)) +
  coord_flip() + labs(title="(a)") +
  xlim(c(1,14)) + theme_bw()+ ylim(c(1,5))+
  geom_point(size=4, col=df_zeta$col) +
  geom_segment(aes(x=x, y=L1, xend=x, yend=U1), size=1,
               col=df_zeta$col) +
  geom_segment(aes(x = x, y = L2, xend = x, yend = U2), size=2,
               col=df_zeta$col) +
  scale_x_continuous(labels=c("M. indica", "M. toxiferum",
                              "R. copallina", "S. terebinthifolia",
                              "T. pubescens", "T. radicans", "T. vernix"),
                     name="Species",
                     breaks=c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5, 13.5))+
  theme(axis.title.x = element_blank())

p2 <- ggplot(df_eta, aes(x=x, y=mean)) +
  coord_flip() + labs(title="(b)") +
  xlim(c(1,14)) + theme_bw()+ ylim(c(-5, 3)) +
  geom_point(size=4, col=df_eta$col) +
  geom_segment(aes(x=x, y=L1, xend=x, yend=U1), size=1,
               col=df_eta$col) +
  geom_segment(aes(x = x, y = L2, xend = x, yend = U2), size=2,
               col=df_eta$col) +
  scale_x_continuous(labels=c("M. indica", "M. toxiferum",
                              "R. copallina", "S. terebinthifolia",
                              "T. pubescens", "T. radicans", "T. vernix"),
                     name="Species",
                     breaks=c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5, 13.5))+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

png("./detection_caterplots.png", width=7,
    height=7, units="in", pointsize=pointsize, res=300
)
grid.arrange(p1,p2, ncol=2, widths=c(1, 0.75))
dev.off()

#Fig 5: Random effects #####

months_df_all <- data.frame(month = c("Jan", "Feb", "Mar", "Apr", "May",
                                      "Jun", "Jul", "Aug", "Sep", "Oct",
                                      "Nov", "Dec"),
                            mean  = rep(NA,12),
                            L1    = rep(NA,12),
                            L2    = rep(NA, 12),
                            U1    = rep(NA,12),
                            U2    = rep(NA,12))

load(paste0(path_to_data, "METTOX/samples_1.RData"))
chain1 <- samples1$samples
rm(samples1)
load(paste0(path_to_data, "METTOX/samples_2.RData"))
chain2 <- samples2$samples
rm(samples2)
mcmc <- mcmc.list(chain1, chain2)
years <- unique(surv_cov$year)
levels(years) #note how they are not in numerical order
collectors <- unique(surv_cov$collector)


mos_all <- MCMCsummary(mcmc, params=c("mos"),
                       probs=c(0.025, 0.16, 0.84, 0.975))[1:12, c(1, 3:6)]
row.names(mos_all) <- c("Jan", "Feb", "Mar", "Apr", "May",
                        "Jun", "Jul", "Aug", "Sep", "Oct",
                        "Nov", "Dec")
colnames(mos_all) <- c("mean", "L1", "L2", "U1", "U2")
mos_all$x <- 1:12


yr_all  <- MCMCsummary(mcmc, params=c("yr"),
                       probs=c(0.025, 0.16, 0.84, 0.975))[, c(1, 3:6)]
row.names(yr_all) <- years
colnames(yr_all) <- c("mean", "L1", "L2", "U1", "U2")
yr_all<-yr_all[order(row.names(yr_all)),] #Sort in ascending order by year
yr_all$x <- 1:length(yr_all$mean)
yr_all$inter <- NA
for(i in 1:length(yr_all$inter)){
  yr_all$inter[i] <- ( yr_all$L1<0) && (yr_all$U2[i]>0)
}
yr_all$col <- mod_cols[3]
for(i in 1:length(yr_all$col)) {
  if(yr_all$inter[i]==1){yr_all$col[i] <- "grey"}
}



obs_all <- MCMCsummary(mcmc, params=c("obs"),
                       probs=c(0.025, 0.16, 0.84, 0.975))[, c(1, 3:6)]
row.names(obs_all) <- collectors
colnames(obs_all) <- c("mean", "L1", "L2", "U1", "U2")
obs_all$x <- 1:length(obs_all$mean)
obs_all$inter <- NA
for(i in 1:length(obs_all$inter)){
  obs_all$inter[i] <- ( obs_all$L1<0) && (obs_all$U2[i]>0)
}
obs_all$col <- rep(mod_cols[3])
for(i in 1:length(obs_all$col)) {
  if(obs_all$inter[i]==1){obs_all$col[i] <- "grey"}
}


#Now for filtered data:
load(paste0(path_to_data, "Metopium toxiferum/surv_cov3.RData"))
surv_cov3 <- subset(surv_cov3, surv_cov3$num_genera > 1 | surv_cov3$h==1)
years_filter <- unique(surv_cov3$year)
collectors_filter <- unique(surv_cov3$collector)
load(paste0(path_to_data, "METTOX/fullmod_filtered_mcmcList.RData"))


mos_filter <- MCMCsummary(mcmc, params=c("mos"),
                          probs=c(0.025, 0.16, 0.84, 0.975))[1:12, c(1, 3:6)]
row.names(mos_filter) <- c("Jan", "Feb", "Mar", "Apr", "May",
                           "Jun", "Jul", "Aug", "Sep", "Oct",
                           "Nov", "Dec")
colnames(mos_filter) <- c("mean", "L1", "L2", "U1", "U2")
mos_filter$x <- 1:length(mos_filter$mean)


yr_filter  <- MCMCsummary(mcmc, params=c("yr"),
                          probs=c(0.025, 0.16, 0.84, 0.975))[, c(1, 3:6)]
row.names(yr_filter) <- years_filter
colnames(yr_filter) <- c("mean", "L1", "L2", "U1", "U2")
yr_filter<-yr_filter[order(row.names(yr_filter)),] #Sort in ascending order by year
yr_filter$x <- 1:length(yr_filter$mean)

yr_filter$inter <- NA
for(i in 1:length(yr_filter$inter)){
  yr_filter$inter[i] <- ( yr_filter$L1<0) && (yr_filter$U2[i]>0)
}
yr_filter$col <- mod_cols[3]
for(i in 1:length(yr_filter$col)) {
  if(yr_filter$inter[i]==1){yr_filter$col[i] <- "grey"}
}


obs_filter <- MCMCsummary(mcmc, params=c("obs"),
                          probs=c(0.025, 0.16, 0.84, 0.975))[, c(1, 3:6)]
row.names(obs_filter) <- collectors_filter
colnames(obs_filter) <- c("mean", "L1", "L2", "U1", "U2")
obs_filter$x <- 1:length(obs_filter$mean)
obs_filter$inter <- NA
for(i in 1:length(obs_filter$inter)){
  obs_filter$inter[i] <- ( obs_filter$L1<0) && (obs_filter$U2[i]>0)
}
obs_filter$col <- rep(mod_cols[4])
for(i in 1:length(obs_filter$col)) {
  if(obs_filter$inter[i]==1){ obs_filter$col[i] <- "grey"}
}


p1 <- ggplot(mos_all, aes(x=x, y=mean)) +
  coord_flip() + labs(title="(a)") +
  theme_bw() + ylim(-10.5,-4) +
  geom_point(size=4, col=mod_cols[3]) +
  geom_segment(aes(x=x, y=L1, xend=x, yend=U2), size=1,
               col=mod_cols_full[1]) +
  geom_segment(aes(x = x, y = L2, xend = x, yend = U1), size=2,
               col=mod_cols_full[1]) +
  scale_x_continuous(labels=c("Jan", "Feb", "Mar", "Apr", "May",
                              "Jun", "Jul", "Aug", "Sep", "Oct",
                              "Nov", "Dec"),
                     name="Month",
                     breaks=1:12)+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12))

p2 <- ggplot(mos_filter, aes(x=x, y=mean)) +
  coord_flip() + labs(title="(b)") +
  theme_bw() + ylim(-10.5,-4) +
  geom_point(size=4, col=mod_cols_full[2]) +
  geom_segment(aes(x=x, y=L1, xend=x, yend=U2), size=1,
               col=mod_cols_full[2]) +
  geom_segment(aes(x = x, y = L2, xend = x, yend = U1), size=2,
               col=mod_cols_full[2]) +
  scale_x_continuous(labels=c("Jan", "Feb", "Mar", "Apr", "May",
                              "Jun", "Jul", "Aug", "Sep", "Oct",
                              "Nov", "Dec"),
                     name="Month",
                     breaks=1:12)+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12))


p3 <-  ggplot(yr_all, aes(x=x, y=mean)) +
  coord_flip() + labs(title="(c)") +
  theme_bw() + ylim(-1.5,1.5) +
  geom_point(size=4, col=yr_all$col) +
  geom_segment(aes(x=x, y=L1, xend=x, yend=U2), size=1,
               col=yr_all$col) +
  geom_segment(aes(x = x, y = L2, xend = x, yend = U1), size=2,
               col=yr_all$col) +
  scale_x_continuous(labels=row.names(yr_all),
                     name="Year",
                     breaks=1:124)+
  theme(axis.title.x = element_blank(),
        plot.title = element_text( size=12),
        axis.title.y = element_text(size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=12))


p4 <-  ggplot(yr_filter, aes(x=x, y=mean)) +
  coord_flip() + labs(title="(e)") +
  theme_bw() + ylim(-1.5,1.5) +
  geom_point(size=4, col=yr_filter$col) +
  geom_segment(aes(x=x, y=L1, xend=x, yend=U2), size=1,
               col=yr_filter$col) +
  geom_segment(aes(x = x, y = L2, xend = x, yend = U1), size=2,
               col=yr_filter$col) +
  scale_x_continuous(labels=row.names(yr_filter),
                     name="Year",
                     breaks=1:length(yr_filter$mean))+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=12))


#Collector Plots
p5 <-  ggplot(obs_all, aes(x=x, y=mean)) +
  coord_flip() + labs(title="(d)") +
  theme_bw() + ylim(-3.5,4.1) +
  geom_point(size=4, col=obs_all$col[obs_all$inter==1], data=subset(obs_all, obs_all$inter==1)) +
  geom_segment(aes(x=x, y=L1, xend=x,
                   yend=U2), size=1,
               col=obs_all$col[obs_all$inter==1],
               data=subset(obs_all, obs_all$inter==1)) +
  geom_segment(aes(x = x[inter==1], y = L2[inter==1], xend = x, yend = U1), size=2,
               col=obs_all$col[obs_all$inter==1],
               data=subset(obs_all, obs_all$inter==1)) +
  geom_point(size=4, col=obs_all$col[obs_all$inter!=1], data=subset(obs_all, obs_all$inter!=1)) +
  geom_segment(aes(x=x, y=L1, xend=x,
                   yend=U2), size=1,
               col=obs_all$col[obs_all$inter!=1],
               data=subset(obs_all, obs_all$inter!=1)) +
  geom_segment(aes(x = x[inter!=1], y = L2[inter!=1], xend = x, yend = U1), size=2,
               col=obs_all$col[obs_all$inter!=1],
               data=subset(obs_all, obs_all$inter!=1)) +
  xlab("Collector") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=12))


star <- row.names(obs_all[obs_all$inter!=1,])
#Judd Patterson, Joseph MDO, F. C. Craigshead
temp <- subset(surv_cov, surv_cov$collector %in% star )
summary(temp)

p6 <-  ggplot(obs_filter, aes(x=x, y=mean)) +
  coord_flip() + labs(title="(f)") +
  theme_bw() + ylim(-3.5,4.1) +
  geom_point(size=4, col=obs_filter$col[obs_filter$inter==1], data=subset(obs_filter, obs_filter$inter==1)) +
  geom_segment(aes(x=x, y=L1, xend=x,
                   yend=U2), size=1,
               col=obs_filter$col[obs_filter$inter==1],
               data=subset(obs_filter, obs_filter$inter==1)) +
  geom_segment(aes(x = x[inter==1], y = L2[inter==1], xend = x, yend = U1), size=2,
               col=obs_filter$col[obs_filter$inter==1],
               data=subset(obs_filter, obs_filter$inter==1)) +
  geom_point(size=4, col=obs_filter$col[obs_filter$inter!=1], data=subset(obs_filter, obs_filter$inter!=1)) +
  geom_segment(aes(x=x, y=L1, xend=x,
                   yend=U2), size=1,
               col=obs_filter$col[obs_filter$inter!=1],
               data=subset(obs_filter, obs_filter$inter!=1)) +
  geom_segment(aes(x = x[inter!=1], y = L2[inter!=1], xend = x, yend = U1), size=2,
               col=obs_filter$col[obs_filter$inter!=1],
               data=subset(obs_filter, obs_filter$inter!=1)) +
  xlab("Collector") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text( size=12),
        axis.title.y = element_text(size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=12))

star2 <- row.names(obs_filter[obs_filter$inter!=1,])
#Ashley M Bradford; D. H. Williams & H.L. Williams"; "G. R. Cooley & R. J. Eaton & J. Monachino
temp2 <- subset(surv_cov, surv_cov$collector %in% star2 )
summary(temp2)

lay <- rbind(c(1,1,2,2),
             c(3,5, 4, 6),
             c(3,5, 4, 6),
             c(3,5, 4, 6))
png("./random_effects.png", width=7,
    height=10, units="in", pointsize=pointsize, res=300
)
grid.arrange(p1, p2, p3, p4, p5, p6, layout_matrix=lay)
dev.off()

#Fig C: Collection bias, other species #####

#MANIND
load(paste0(path_to_data, "Mangifera indica/surv_cov.RData"))
surv_cov$prev_detected <- factor(surv_cov$prev_detected)
levels(surv_cov$prev_detected) <- c("First Encounter", "Repeat Encounter")
temp <- subset(deduplicated, deduplicated$species == "Mangifera indica")
counts <- as.data.frame(table(temp$month))
decades <- as.data.frame(table(temp$decade))
citSci <- as.data.frame(table(temp$citizenScience))
detStat <- as.data.frame(table(surv_cov$prev_detected))

p1 <- ggplot(counts, aes(x=Var1, y= Freq)) +
  geom_bar(stat="identity") +
  ggtitle("(a) ") +
  xlab("Month") + ylab("Number of Records") +
  theme(text = element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90))

p2 <- ggplot(decades, aes(x=Var1, y= Freq)) +
  geom_bar(stat="identity") +
  ggtitle("(b) ") +
  xlab("Decade") + ylab("Number of Records") +
  theme(text = element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90))

p3 <- ggplot(citSci, aes(x=Var1, y= Freq)) +
  geom_bar(stat="identity") +
  ggtitle("(c) ") +
  xlab("Collection Type") + ylab("Number of Records") +
  theme(text = element_text(size=10))

p4 <- ggplot(detStat, aes(x=Var1, y= Freq)) +
  geom_bar(stat="identity") +
  ggtitle("(d) ") +
  xlab("Collection Type") + ylab("Number of Records") +
  theme(text = element_text(size=10))

png("./collection_bias_MANIND.png", width=8,
    height=10, units="in", pointsize=pointsize, res=300
)

grid.arrange(p1,p2,p3,p4, ncol=2)
dev.off()

#2. RHUCOP
load(paste0(path_to_data, "Rhus copallina/surv_cov.RData"))
surv_cov$prev_detected <- factor(surv_cov$prev_detected)
levels(surv_cov$prev_detected) <- c("First Encounter", "Repeat Encounter")
temp <- subset(deduplicated, deduplicated$species == "Rhus copallina")
counts <- as.data.frame(table(temp$month))
decades <- as.data.frame(table(temp$decade))
citSci <- as.data.frame(table(temp$citizenScience))
detStat <- as.data.frame(table(surv_cov$prev_detected))

p1 <- ggplot(counts, aes(x=Var1, y= Freq)) +
  geom_bar(stat="identity") +
  ggtitle("(a) ") +
  xlab("Month") + ylab("Number of Records") +
  theme(text = element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90))

p2 <- ggplot(decades, aes(x=Var1, y= Freq)) +
  geom_bar(stat="identity") +
  ggtitle("(b) ") +
  xlab("Decade") + ylab("Number of Records") +
  theme(text = element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90))

p3 <- ggplot(citSci, aes(x=Var1, y= Freq)) +
  geom_bar(stat="identity") +
  ggtitle("(c) ") +
  xlab("Collection Type") + ylab("Number of Records") +
  theme(text = element_text(size=10))

p4 <- ggplot(detStat, aes(x=Var1, y= Freq)) +
  geom_bar(stat="identity") +
  ggtitle("(d) ") +
  xlab("Collection Type") + ylab("Number of Records") +
  theme(text = element_text(size=10))

png("./collection_bias_RHUCOP.png", width=8,
    height=10, units="in", pointsize=pointsize, res=300
)

grid.arrange(p1,p2,p3,p4, ncol=2)
dev.off()

#3. SCHTER
load(paste0(path_to_data, "Schinus terebinthifolia/surv_cov.RData"))
surv_cov$prev_detected <- factor(surv_cov$prev_detected)
levels(surv_cov$prev_detected) <- c("First Encounter", "Repeat Encounter")
temp <- subset(deduplicated, deduplicated$species == "Schinus terebinthifolia")
counts <- as.data.frame(table(temp$month))
decades <- as.data.frame(table(temp$decade))
citSci <- as.data.frame(table(temp$citizenScience))
detStat <- as.data.frame(table(surv_cov$prev_detected))

p1 <- ggplot(counts, aes(x=Var1, y= Freq)) +
  geom_bar(stat="identity") +
  ggtitle("(a) ") +
  xlab("Month") + ylab("Number of Records") +
  theme(text = element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90))

p2 <- ggplot(decades, aes(x=Var1, y= Freq)) +
  geom_bar(stat="identity") +
  ggtitle("(b) ") +
  xlab("Decade") + ylab("Number of Records") +
  theme(text = element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90))

p3 <- ggplot(citSci, aes(x=Var1, y= Freq)) +
  geom_bar(stat="identity") +
  ggtitle("(c) ") +
  xlab("Collection Type") + ylab("Number of Records") +
  theme(text = element_text(size=10))

p4 <- ggplot(detStat, aes(x=Var1, y= Freq)) +
  geom_bar(stat="identity") +
  ggtitle("(d) ") +
  xlab("Collection Type") + ylab("Number of Records") +
  theme(text = element_text(size=10))

png("./collection_bias_SCHTER.png", width=8,
    height=10, units="in", pointsize=pointsize, res=300
)

grid.arrange(p1,p2,p3,p4, ncol=2)
dev.off()

#4. TOXPUB
load(paste0(path_to_data, "Toxicodendron pubescens/surv_cov.RData"))
surv_cov$prev_detected <- factor(surv_cov$prev_detected)
levels(surv_cov$prev_detected) <- c("First Encounter", "Repeat Encounter")
temp <- subset(deduplicated, deduplicated$species == "Toxicodendron pubescens")
counts <- as.data.frame(table(temp$month))
decades <- as.data.frame(table(temp$decade))
citSci <- as.data.frame(table(temp$citizenScience))
detStat <- as.data.frame(table(surv_cov$prev_detected))

p1 <- ggplot(counts, aes(x=Var1, y= Freq)) +
  geom_bar(stat="identity") +
  ggtitle("(a) ") +
  xlab("Month") + ylab("Number of Records") +
  theme(text = element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90))

p2 <- ggplot(decades, aes(x=Var1, y= Freq)) +
  geom_bar(stat="identity") +
  ggtitle("(b) ") +
  xlab("Decade") + ylab("Number of Records") +
  theme(text = element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90))

p3 <- ggplot(citSci, aes(x=Var1, y= Freq)) +
  geom_bar(stat="identity") +
  ggtitle("(c) ") +
  xlab("Collection Type") + ylab("Number of Records") +
  theme(text = element_text(size=10))

p4 <- ggplot(detStat, aes(x=Var1, y= Freq)) +
  geom_bar(stat="identity") +
  ggtitle("(d) ") +
  xlab("Collection Type") + ylab("Number of Records") +
  theme(text = element_text(size=10))

png("./collection_bias_TOXPUB.png", width=8,
    height=10, units="in", pointsize=pointsize, res=300
)

grid.arrange(p1,p2,p3,p4, ncol=2)
dev.off()

#5. TOXRAD
load(paste0(path_to_data, "Toxicodendron radicans/surv_cov.RData"))
surv_cov$prev_detected <- factor(surv_cov$prev_detected)
levels(surv_cov$prev_detected) <- c("First Encounter", "Repeat Encounter")
temp <- subset(deduplicated, deduplicated$species == "Toxicodendron radicans")
counts <- as.data.frame(table(temp$month))
decades <- as.data.frame(table(temp$decade))
citSci <- as.data.frame(table(temp$citizenScience))
detStat <- as.data.frame(table(surv_cov$prev_detected))

p1 <- ggplot(counts, aes(x=Var1, y= Freq)) +
  geom_bar(stat="identity") +
  ggtitle("(a) ") +
  xlab("Month") + ylab("Number of Records") +
  theme(text = element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90))

p2 <- ggplot(decades, aes(x=Var1, y= Freq)) +
  geom_bar(stat="identity") +
  ggtitle("(b) ") +
  xlab("Decade") + ylab("Number of Records") +
  theme(text = element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90))

p3 <- ggplot(citSci, aes(x=Var1, y= Freq)) +
  geom_bar(stat="identity") +
  ggtitle("(c) ") +
  xlab("Collection Type") + ylab("Number of Records") +
  theme(text = element_text(size=10))

p4 <- ggplot(detStat, aes(x=Var1, y= Freq)) +
  geom_bar(stat="identity") +
  ggtitle("(d) ") +
  xlab("Collection Type") + ylab("Number of Records") +
  theme(text = element_text(size=10))

png("./collection_bias_TOXRAD.png", width=8,
    height=10, units="in", pointsize=pointsize, res=300
)

grid.arrange(p1,p2,p3,p4, ncol=2)
dev.off()

#6. TOXVER
load(paste0(path_to_data, "Toxicodendron vernix/surv_cov.RData"))
surv_cov$prev_detected <- factor(surv_cov$prev_detected)
levels(surv_cov$prev_detected) <- c("First Encounter", "Repeat Encounter")
temp <- subset(deduplicated, deduplicated$species == "Toxicodendron vernix")
counts <- as.data.frame(table(temp$month))
decades <- as.data.frame(table(temp$decade))
citSci <- as.data.frame(table(temp$citizenScience))
detStat <- as.data.frame(table(surv_cov$prev_detected))

p1 <- ggplot(counts, aes(x=Var1, y= Freq)) +
  geom_bar(stat="identity") +
  ggtitle("(a) ") +
  xlab("Month") + ylab("Number of Records") +
  theme(text = element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90))

p2 <- ggplot(decades, aes(x=Var1, y= Freq)) +
  geom_bar(stat="identity") +
  ggtitle("(b) ") +
  xlab("Decade") + ylab("Number of Records") +
  theme(text = element_text(size=10)) +
  theme(axis.text.x = element_text(angle=90))

p3 <- ggplot(citSci, aes(x=Var1, y= Freq)) +
  geom_bar(stat="identity") +
  ggtitle("(c) ") +
  xlab("Collection Type") + ylab("Number of Records") +
  theme(text = element_text(size=10))

p4 <- ggplot(detStat, aes(x=Var1, y= Freq)) +
  geom_bar(stat="identity") +
  ggtitle("(d) ") +
  xlab("Collection Type") + ylab("Number of Records") +
  theme(text = element_text(size=10))

png("./collection_bias_TOXVER.png", width=8,
    height=10, units="in", pointsize=pointsize, res=300
)

grid.arrange(p1,p2,p3,p4, ncol=2)
dev.off()

#Fig D: Occupancy Caterplots #####
#Posterior caterplots of occupancy covariates
#beta: human population density
#gamma: minimum temperature
#delta: county area

speciesList <- c("MANIND", "METTOX", "RHUCOP", "SCHTER", "TOXPUB", "TOXRAD", "TOXVER")

df_beta <- data.frame(species = rep(c("MANIND", "METTOX", "RHUCOP", "SCHTER", "TOXPUB", "TOXRAD", "TOXVER"), each= 2),
                      model = rep(c("full", "filter"), 7),
                      x = 1:14,
                      mean = rep(NA,14),
                      L1 = rep(NA, 14),
                      L2 = rep(NA, 14),
                      U2 = rep(NA,14),
                      U1 = rep(NA,14))

df_gamma <- data.frame(species = rep(c("MANIND", "METTOX", "RHUCOP", "SCHTER", "TOXPUB", "TOXRAD", "TOXVER"), each= 2),
                       model = rep(c("full", "filter"), 7),
                       x = 1:14,
                       mean = rep(NA,14),
                       L1 = rep(NA, 14),
                       L2 = rep(NA, 14),
                       U2 = rep(NA,14),
                       U1 = rep(NA,14))
df_delta <- data.frame(species = rep(c("MANIND", "METTOX", "RHUCOP", "SCHTER", "TOXPUB", "TOXRAD", "TOXVER"), each= 2),
                       model = rep(c("full", "filter"), 7),
                       x = 1:14,
                       mean = rep(NA,14),
                       L1 = rep(NA, 14),
                       L2 = rep(NA, 14),
                       U2 = rep(NA,14),
                       U1 = rep(NA,14))



pos <- c(1, 3, 5, 7, 9, 11, 13)


for (i in 2:7) {
  load(paste0(path_to_data, speciesList[i], "/samples_1.RData"))
  chain1 <- samples1$samples
  rm(samples1)
  load(paste0(path_to_data, speciesList[i], "/samples_2.RData"))
  chain2 <- samples2$samples
  rm(samples2)
  mcmc <- mcmc.list(chain1, chain2)


  thing_all_beta <- MCMCsummary(mcmc, params=c("beta"),
                                probs=c(0.025, 0.16, 0.84, 0.975))[1, c(1, 3:6)]
  df_beta[pos[i], 4:8] <- thing_all_beta
  rm(thing_all_beta)

  thing_all_gamma <- MCMCsummary(mcmc, params=c("gamma"),
                                 probs=c(0.025, 0.16, 0.84, 0.975))[1, c(1, 3:6)]
  df_gamma[pos[i], 4:8] <- thing_all_gamma
  rm(thing_all_gamma)

  thing_all_delta <- MCMCsummary(mcmc, params=c("delta"),
                                 probs=c(0.025, 0.16, 0.84, 0.975))[1, c(1, 3:6)]
  df_delta[pos[i], 4:8] <- thing_all_delta
  rm(thing_all_delta)

  rm(mcmc)

  load(paste0(path_to_data, speciesList[i], "/fullmod_filtered_mcmcList.RData"))

  thing_filter_beta <- MCMCsummary(mcmc, params=c("beta"),
                                   probs=c(0.025, 0.16, 0.84, 0.975))[1, c(1, 3:6)]
  df_beta[(pos[i]+1), 4:8] <- thing_filter_beta
  rm(thing_filter_beta)

  thing_filter_gamma <- MCMCsummary(mcmc, params=c("gamma"),
                                    probs=c(0.025, 0.16, 0.84, 0.975))[1, c(1, 3:6)]
  df_gamma[(pos[i]+1), 4:8] <- thing_filter_gamma
  rm(thing_filter_gamma)

  thing_filter_delta <- MCMCsummary(mcmc, params=c("delta"),
                                    probs=c(0.025, 0.16, 0.84, 0.975))[1, c(1, 3:6)]
  df_delta[(pos[i]+1), 4:8] <- thing_filter_delta
  rm(thing_filter_delta)
  rm(mcmc)
}


min(df_beta$L1) # -4.51
max(df_beta$U1) # 4.67
min(df_gamma$L1) # -4.98
max(df_gamma$U1) # 4.97
min(df_delta$L1) # -2.72
max(df_delta$U1) # 4.82

#Adjust color scheme so that CI that overlap zero are greyed out
#beta
df_beta$inter <- NA
for(i in 1:length(df_beta$inter)){
  df_beta$inter[i] <- ( df_beta$L1[i]<0) && (df_beta$U2[i]>0)
}
df_beta$col <- c(mod_cols[3], mod_cols[4])
for(i in 1:length(df_beta$col)) {
  if(df_beta$inter[i]==1){df_beta$col[i] <- "grey"}
}

#gamma
df_gamma$inter <- NA
for(i in 1:length(df_gamma$inter)){
  df_gamma$inter[i] <- ( df_gamma$L1[i]<0) && (df_gamma$U2[i]>0)
}
df_gamma$col <- c(mod_cols[3], mod_cols[4])
for(i in 1:length(df_gamma$col)) {
  if(df_gamma$inter[i]==1){df_gamma$col[i] <- "grey"}
}

#delta
df_delta$inter <- NA
for(i in 1:length(df_delta$inter)){
  df_delta$inter[i] <- ( df_delta$L1[i]<0) && (df_delta$U2[i]>0)
}
df_delta$col <- c(mod_cols[3], mod_cols[4])
for(i in 1:length(df_delta$col)) {
  if(df_delta$inter[i]==1){df_delta$col[i] <- "grey"}
}



save(df_beta, file=paste0(path_to_data, "crossval_beta.RData"))
save(df_gamma, file=paste0(path_to_data, "crossval_gamma.RData"))
save(df_delta, file=paste0(path_to_data, "crossval_delta.RData"))

p1 <- ggplot(df_beta, aes(x=x, y=mean)) +
  coord_flip() + labs(title="(a)") +
  xlim(c(1,14)) + theme_bw()+ ylim(c(-5,5))+
  geom_point(size=4, col=df_beta$col) +
  geom_segment(aes(x=x, y=L1, xend=x, yend=U1), size=1,
               col=df_beta$col) +
  geom_segment(aes(x = x, y = L2, xend = x, yend = U2), size=2,
               col=df_beta$col) +
  scale_x_continuous(labels=c("M. indica", "M. toxiferum",
                              "R. copallina", "S. terebinthifolia",
                              "T. pubescens", "T. radicans", "T. vernix"),
                     name="Species",
                     breaks=c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5, 13.5))+
  theme(axis.title.x = element_blank())

p2 <- ggplot(df_gamma, aes(x=x, y=mean)) +
  coord_flip() + labs(title="(b)") +
  xlim(c(1,14)) + theme_bw()+ ylim(c(-5, 5)) +
  geom_point(size=4, col=df_gamma$col) +
  geom_segment(aes(x=x, y=L1, xend=x, yend=U1), size=1,
               col=df_gamma$col) +
  geom_segment(aes(x = x, y = L2, xend = x, yend = U2), size=2,
               col=df_gamma$col) +
  scale_x_continuous(labels=c("M. indica", "M. toxiferum",
                              "R. copallina", "S. terebinthifolia",
                              "T. pubescens", "T. radicans", "T. vernix"),
                     name="Species",
                     breaks=c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5, 13.5))+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

p3 <- ggplot(df_delta, aes(x=x, y=mean)) +
  coord_flip() + labs(title="(c)") +
  xlim(c(1,14)) + theme_bw()+ ylim(c(-3, 5)) +
  geom_point(size=4, col=df_delta$col) +
  geom_segment(aes(x=x, y=L1, xend=x, yend=U1), size=1,
               col=df_delta$col) +
  geom_segment(aes(x = x, y = L2, xend = x, yend = U2), size=2,
               col=df_delta$col) +
  scale_x_continuous(labels=c("M. indica", "M. toxiferum",
                              "R. copallina", "S. terebinthifolia",
                              "T. pubescens", "T. radicans", "T. vernix"),
                     name="Species",
                     breaks=c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5, 13.5))+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

png("./occupancy_caterplots.png", width=7,
    height=7, units="in", pointsize=pointsize, res=300
)
grid.arrange(p1,p2, p3, ncol=3, widths=c(1, 0.67, 0.67))
dev.off()

#Fig E: Random effects for remaining species #####
speciesCodeList <- c("MANIND", "RHUCOP", "SCHTER", "TOXPUB", "TOXRAD", "TOXVER")
speciesNameList <- c("Mangifera indica",
                     "Rhus copallina",
                     "Schinus terebinthifolia",
                     "Toxicodendron pubescens",
                     "Toxicodendron radicans",
                     "Toxicodendron vernix")

#Have to define axis limits for species
#mos_lim
min(mos_all$L1, mos_filter$L1)
max(mos_all$U2, mos_filter$U2)
#MANIND c(-27, -1)
mos_lim <- list(c(-27, -1), #MANIND
                c(-9,-4),  #RHUCOP
                c(-8, -2), #SCHTER
                c(-15, 1), #TOXPUB
                c(-11, -5), #TOXRAD
                c(-9,1)) #TOXVER

#yr_lim
min(yr_all$L1, yr_filter$L1)
max(yr_all$U2, yr_filter$U2)
#MANIND c(-2, 2.5)
yr_lim <- list(c(-2, 2.5), #MANIND
               c(-1,1 ), #RHUCOP
               c(-1.5, 1.5), #SCHTER
               c(-1, 2), #TOXPUB
               c(-2, 2), #TOXRAD
               c(-1, 2))

#obs_lim
min(obs_all$L1, obs_filter$L1)
max(obs_all$U2, obs_filter$U2)
#MANIND c(-5, 10)
obs_lim <- list(c(-5,10), #MANIND
                c(-4,4 ), #RHUCOP
                c(-4.5, 5), #SCHTER
                c(-3,3), #TOXPUB
                c(-5, 6), #TOXRAD
                c(-4, 5))

for (s in 1:length(speciesCodeList)) {



  load(paste0(path_to_data, speciesCodeList[s], "/samples_1.RData"))
  chain1 <- samples1$samples
  rm(samples1)
  load(paste0(path_to_data, speciesCodeList[s], "/samples_2.RData"))
  chain2 <- samples2$samples
  rm(samples2)
  mcmc <- mcmc.list(chain1, chain2)
  years <- unique(surv_cov$year)
  levels(years) #note how they are not in numerical order
  collectors <- unique(surv_cov$collector)


  mos_all <- MCMCsummary(mcmc, params=c("mos"),
                         probs=c(0.025, 0.16, 0.84, 0.975))[1:12, c(1, 3:6)]
  row.names(mos_all) <- c("Jan", "Feb", "Mar", "Apr", "May",
                          "Jun", "Jul", "Aug", "Sep", "Oct",
                          "Nov", "Dec")
  colnames(mos_all) <- c("mean", "L1", "L2", "U1", "U2")
  mos_all$x <- 1:12


  yr_all  <- MCMCsummary(mcmc, params=c("yr"),
                         probs=c(0.025, 0.16, 0.84, 0.975))[, c(1, 3:6)]
  row.names(yr_all) <- years
  colnames(yr_all) <- c("mean", "L1", "L2", "U1", "U2")
  yr_all<-yr_all[order(row.names(yr_all)),] #Sort in ascending order by year
  yr_all$x <- 1:length(yr_all$mean)
  yr_all$inter <- NA
  for(i in 1:length(yr_all$inter)){
    yr_all$inter[i] <- ( yr_all$L1<0) && (yr_all$U2[i]>0)
  }
  yr_all$col <- mod_cols[3]
  for(i in 1:length(yr_all$col)) {
    if(yr_all$inter[i]==1){yr_all$col[i] <- "grey"}
  }



  obs_all <- MCMCsummary(mcmc, params=c("obs"),
                         probs=c(0.025, 0.16, 0.84, 0.975))[, c(1, 3:6)]
  row.names(obs_all) <- collectors
  colnames(obs_all) <- c("mean", "L1", "L2", "U1", "U2")
  obs_all$x <- 1:length(obs_all$mean)
  obs_all$inter <- NA
  for(i in 1:length(obs_all$inter)){
    obs_all$inter[i] <- ( obs_all$L1<0) && (obs_all$U2[i]>0)
  }
  obs_all$col <- rep(mod_cols[3])
  for(i in 1:length(obs_all$col)) {
    if(obs_all$inter[i]==1){obs_all$col[i] <- "grey"}
  }


  #Now for filtered data:
  load(paste0(path_to_data, speciesNameList[s], "/surv_cov3.RData"))
  surv_cov3 <- subset(surv_cov3, surv_cov3$num_genera > 1 | surv_cov3$h==1)
  years_filter <- unique(surv_cov3$year)
  collectors_filter <- unique(surv_cov3$collector)
  load(paste0(path_to_data, speciesCodeList[s], "/fullmod_filtered_mcmcList.RData"))


  mos_filter <- MCMCsummary(mcmc, params=c("mos"),
                            probs=c(0.025, 0.16, 0.84, 0.975))[1:12, c(1, 3:6)]
  row.names(mos_filter) <- c("Jan", "Feb", "Mar", "Apr", "May",
                             "Jun", "Jul", "Aug", "Sep", "Oct",
                             "Nov", "Dec")
  colnames(mos_filter) <- c("mean", "L1", "L2", "U1", "U2")
  mos_filter$x <- 1:length(mos_filter$mean)


  yr_filter  <- MCMCsummary(mcmc, params=c("yr"),
                            probs=c(0.025, 0.16, 0.84, 0.975))[, c(1, 3:6)]
  row.names(yr_filter) <- years_filter
  colnames(yr_filter) <- c("mean", "L1", "L2", "U1", "U2")
  yr_filter<-yr_filter[order(row.names(yr_filter)),] #Sort in ascending order by year
  yr_filter$x <- 1:length(yr_filter$mean)

  yr_filter$inter <- NA
  for(i in 1:length(yr_filter$inter)){
    yr_filter$inter[i] <- ( yr_filter$L1<0) && (yr_filter$U2[i]>0)
  }
  yr_filter$col <- mod_cols[3]
  for(i in 1:length(yr_filter$col)) {
    if(yr_filter$inter[i]==1){yr_filter$col[i] <- "grey"}
  }


  obs_filter <- MCMCsummary(mcmc, params=c("obs"),
                            probs=c(0.025, 0.16, 0.84, 0.975))[, c(1, 3:6)]
  row.names(obs_filter) <- collectors_filter
  colnames(obs_filter) <- c("mean", "L1", "L2", "U1", "U2")
  obs_filter$x <- 1:length(obs_filter$mean)
  obs_filter$inter <- NA
  for(i in 1:length(obs_filter$inter)){
    obs_filter$inter[i] <- ( obs_filter$L1<0) && (obs_filter$U2[i]>0)
  }
  obs_filter$col <- rep(mod_cols[4])
  for(i in 1:length(obs_filter$col)) {
    if(obs_filter$inter[i]==1){ obs_filter$col[i] <- "grey"}
  }


  save(obs_all, file=paste0(path_to_data, speciesNameList[s], "/", speciesCodeList[s], "_obs_all.RData"))
  save(obs_filter, file=paste0(path_to_data, speciesNameList[s], "/", speciesCodeList[s], "_obs_filter.RData"))

  save(mos_all, file=paste0(path_to_data, speciesNameList[s], "/", speciesCodeList[s], "_mos_all.RData"))
  save(mos_filter, file=paste0(path_to_data, speciesNameList[s], "/", speciesCodeList[s], "_mos_filter.RData"))

  save(yr_all, file=paste0(path_to_data, speciesNameList[s], "/", speciesCodeList[s], "_yr_all.RData"))
  save(yr_filter, file=paste0(path_to_data, speciesNameList[s], "/", speciesCodeList[s], "_yr_filter.RData"))

}



p1 <- ggplot(mos_all, aes(x=x, y=mean)) +
  coord_flip() + labs(title="(a)") +
  theme_bw() + ylim(mos_lim[[s]]) +
  geom_point(size=4, col=mod_cols[3]) +
  geom_segment(aes(x=x, y=L1, xend=x, yend=U2), size=1,
               col=mod_cols_full[1]) +
  geom_segment(aes(x = x, y = L2, xend = x, yend = U1), size=2,
               col=mod_cols_full[1]) +
  scale_x_continuous(labels=c("Jan", "Feb", "Mar", "Apr", "May",
                              "Jun", "Jul", "Aug", "Sep", "Oct",
                              "Nov", "Dec"),
                     name="Month",
                     breaks=1:12)+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12))

p2 <- ggplot(mos_filter, aes(x=x, y=mean)) +
  coord_flip() + labs(title="(b)") +
  theme_bw() + ylim(mos_lim[[s]]) +
  geom_point(size=4, col=mod_cols_full[2]) +
  geom_segment(aes(x=x, y=L1, xend=x, yend=U2), size=1,
               col=mod_cols_full[2]) +
  geom_segment(aes(x = x, y = L2, xend = x, yend = U1), size=2,
               col=mod_cols_full[2]) +
  scale_x_continuous(labels=c("Jan", "Feb", "Mar", "Apr", "May",
                              "Jun", "Jul", "Aug", "Sep", "Oct",
                              "Nov", "Dec"),
                     name="Month",
                     breaks=1:12)+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12))


p3 <-  ggplot(yr_all, aes(x=x, y=mean)) +
  coord_flip() + labs(title="(c)") +
  theme_bw() + ylim(yr_lim[[s]]) +
  geom_point(size=4, col=yr_all$col) +
  geom_segment(aes(x=x, y=L1, xend=x, yend=U2), size=1,
               col=yr_all$col) +
  geom_segment(aes(x = x, y = L2, xend = x, yend = U1), size=2,
               col=yr_all$col) +
  scale_x_continuous(labels=row.names(yr_all),
                     name="Year",
                     breaks=1:124)+
  theme(axis.title.x = element_blank(),
        plot.title = element_text( size=12),
        axis.title.y = element_text(size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=12))


p4 <-  ggplot(yr_filter, aes(x=x, y=mean)) +
  coord_flip() + labs(title="(e)") +
  theme_bw() + ylim(yr_lim[[s]]) +
  geom_point(size=4, col=yr_filter$col) +
  geom_segment(aes(x=x, y=L1, xend=x, yend=U2), size=1,
               col=yr_filter$col) +
  geom_segment(aes(x = x, y = L2, xend = x, yend = U1), size=2,
               col=yr_filter$col) +
  scale_x_continuous(labels=row.names(yr_filter),
                     name="Year",
                     breaks=1:length(yr_filter$mean))+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=12))


#Collector Plots
p5 <-  ggplot(obs_all, aes(x=x, y=mean)) +
  coord_flip() + labs(title="(d)") +
  theme_bw() + ylim(obs_lim[[s]]) +
  geom_point(size=4, col=obs_all$col[obs_all$inter==1], data=subset(obs_all, obs_all$inter==1)) +
  geom_segment(aes(x=x, y=L1, xend=x,
                   yend=U2), size=1,
               col=obs_all$col[obs_all$inter==1],
               data=subset(obs_all, obs_all$inter==1)) +
  geom_segment(aes(x = x[inter==1], y = L2[inter==1], xend = x, yend = U1), size=2,
               col=obs_all$col[obs_all$inter==1],
               data=subset(obs_all, obs_all$inter==1)) +
  geom_point(size=4, col=obs_all$col[obs_all$inter!=1], data=subset(obs_all, obs_all$inter!=1)) +
  geom_segment(aes(x=x, y=L1, xend=x,
                   yend=U2), size=1,
               col=obs_all$col[obs_all$inter!=1],
               data=subset(obs_all, obs_all$inter!=1)) +
  geom_segment(aes(x = x[inter!=1], y = L2[inter!=1], xend = x, yend = U1), size=2,
               col=obs_all$col[obs_all$inter!=1],
               data=subset(obs_all, obs_all$inter!=1)) +
  xlab("Collector") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=12))


p6 <-  ggplot(obs_filter, aes(x=x, y=mean)) +
  coord_flip() + labs(title="(f)") +
  theme_bw() + ylim(obs_lim[[s]]) +
  geom_point(size=4, col=obs_filter$col[obs_filter$inter==1], data=subset(obs_filter, obs_filter$inter==1)) +
  geom_segment(aes(x=x, y=L1, xend=x,
                   yend=U2), size=1,
               col=obs_filter$col[obs_filter$inter==1],
               data=subset(obs_filter, obs_filter$inter==1)) +
  geom_segment(aes(x = x[inter==1], y = L2[inter==1], xend = x, yend = U1), size=2,
               col=obs_filter$col[obs_filter$inter==1],
               data=subset(obs_filter, obs_filter$inter==1)) +
  geom_point(size=4, col=obs_filter$col[obs_filter$inter!=1], data=subset(obs_filter, obs_filter$inter!=1)) +
  geom_segment(aes(x=x, y=L1, xend=x,
                   yend=U2), size=1,
               col=obs_filter$col[obs_filter$inter!=1],
               data=subset(obs_filter, obs_filter$inter!=1)) +
  geom_segment(aes(x = x[inter!=1], y = L2[inter!=1], xend = x, yend = U1), size=2,
               col=obs_filter$col[obs_filter$inter!=1],
               data=subset(obs_filter, obs_filter$inter!=1)) +
  xlab("Collector") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text( size=12),
        axis.title.y = element_text(size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=12))


lay <- rbind(c(1,1,2,2),
             c(3,5, 4, 6),
             c(3,5, 4, 6),
             c(3,5, 4, 6))
png(paste0("./", speciesCodeList[s], "_random_effects.png"), width=7,
    height=10, units="in", pointsize=pointsize, res=300
)
grid.arrange(p1, p2, p3, p4, p5, p6, layout_matrix=lay)
dev.off()

}
