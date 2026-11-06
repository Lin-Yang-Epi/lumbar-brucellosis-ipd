# =============================================================================
# 03_tableS1_fig2_map.R
# Table S1. Countries of the included patients and their endemic classification
# Fig. 2.  World map: patients per country, countries shaded by current endemic status
#
# Run from the project root (the folder that contains data/ and scripts/).
# Input : data/clean/ipd_clean.csv
#         data/raw/endemic_classification.csv
# Output: outputs/tables/TableS1_countries.csv
#         outputs/figures/Fig2_world_map.pdf / .png / .tiff
# Packages: ggplot2, sf, rnaturalearth, rnaturalearthdata, patchwork
# =============================================================================

library(ggplot2)
library(sf)
library(rnaturalearth)

ipd <- read.csv(file.path("data", "clean", "ipd_clean.csv"),
                na.strings = "NA", stringsAsFactors = FALSE)
ec  <- read.csv(file.path("data", "raw", "endemic_classification.csv"),
                stringsAsFactors = FALSE, check.names = FALSE)
ec$Region <- trimws(ec$Region)

# -----------------------------------------------------------------------------
# 1. Country lookup: map code (Natural Earth ADM0_A3) and bubble position
#    French Guiana is drawn as part of France on the map; its bubble sits on
#    French Guiana itself. Russia's bubble sits on the Republic of Dagestan, where the
#    reported patient (Burgasova 2026) was diagnosed; the whole country is shaded.
# -----------------------------------------------------------------------------
lookup <- read.csv(text = "
Region,adm0_a3,lon,lat
China,CHN,104.0,35.0
Turkey,TUR,35.0,39.0
Spain,ESP,-3.7,40.2
USA,USA,-98.0,39.0
India,IND,79.0,22.0
Greece,GRC,22.0,39.5
Iran,IRN,53.0,32.5
Saudi Arabia,SAU,45.0,24.0
UK,GBR,-2.0,53.5
Germany,DEU,10.4,51.0
Korea,KOR,127.8,36.3
Italy,ITA,12.6,42.8
Australia,AUS,134.0,-25.5
Bosnia and Herzegovina,BIH,17.8,44.2
Canada,CAN,-106.0,56.0
Kazakhstan,KAZ,67.0,48.0
Malaysia,MYS,102.3,4.0
Peru,PER,-75.0,-9.5
Tunisia,TUN,9.5,34.0
Algeria,DZA,2.6,28.0
Brazil,BRA,-52.0,-10.0
Czech Republic,CZE,15.5,49.8
France,FRA,2.4,46.6
French Guiana,FRA,-53.1,4.0
Israel,ISR,34.9,31.4
Japan,JPN,138.0,36.5
Jordan,JOR,36.5,31.0
Kuwait,KWT,47.6,29.3
Lebanon,LBN,35.9,33.9
Morocco,MAR,-6.5,31.8
Netherlands,NLD,5.3,52.2
New Zealand,NZL,172.5,-41.5
Norway,NOR,9.0,61.0
Oman,OMN,57.0,21.0
Portugal,PRT,-8.2,39.6
Qatar,QAT,51.2,25.3
Somalia,SOM,46.0,6.0
South Africa,ZAF,24.0,-29.0
Sri Lanka,LKA,80.7,7.9
Syria,SYR,38.5,35.0
Thailand,THA,101.0,15.0
Armenia,ARM,45.0,40.2
Colombia,COL,-74.0,4.5
Russia,RUS,47.0,42.6
", stringsAsFactors = FALSE, strip.white = TRUE)

missing_lookup <- setdiff(unique(ipd$Region), lookup$Region)
if (length(missing_lookup)) stop("Add these countries to the lookup table: ",
                                 paste(missing_lookup, collapse = ", "))
missing_ec <- setdiff(unique(ipd$Region), ec$Region)
if (length(missing_ec)) stop("Not in endemic_classification.csv: ",
                             paste(missing_ec, collapse = ", "))

# -----------------------------------------------------------------------------
# 2. Table S1
# -----------------------------------------------------------------------------
counts <- as.data.frame(table(Region = ipd$Region), stringsAsFactors = FALSE)
names(counts)[2] <- "n"
tab <- merge(counts, ec, by = "Region", all.x = TRUE)
if ("N_patients" %in% names(tab)) {
  off <- tab[!is.na(tab$N_patients) & tab$N_patients != tab$n, c("Region", "n", "N_patients")]
  if (nrow(off)) warning("N_patients in endemic_classification.csv differs from the data:\n",
                         paste(capture.output(print(off, row.names = FALSE)), collapse = "\n"))
}
tab$pct <- round(100 * tab$n / nrow(ipd), 1)
tab <- tab[order(-tab$n, tab$Region), ]

table_s1 <- data.frame(
  Country                         = tab$Region,
  `Patients, n`                   = tab$n,
  `Patients, %`                   = tab$pct,
  `Current status`                = tab$Current_class_proposed,
  `Evidence level`                = tab$Current_evidence,
  `Source or note`                = tab$Current_source_or_note,
  `Pappas 2006 incidence (per million)` = tab$Pappas2006_incidence_per_million,
  `Pappas 2006 status`            = tab$Pappas2006_class,
  `Status differs`                = ifelse(tab$Differs %in% "Y", "Yes", ""),
  check.names = FALSE, stringsAsFactors = FALSE)
dir.create(file.path("outputs", "tables"), recursive = TRUE, showWarnings = FALSE)
write.csv(table_s1, file.path("outputs", "tables", "TableS1_countries.csv"),
          row.names = FALSE, na = "")

cat("Countries:", nrow(tab), "| patients:", sum(tab$n), "\n")
by_status <- aggregate(n ~ Current_class_proposed, data = tab, FUN = sum)
by_status$countries <- as.vector(table(tab$Current_class_proposed)[by_status$Current_class_proposed])
print(by_status, row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. Map data
# -----------------------------------------------------------------------------
world <- ne_countries(scale = "medium", returnclass = "sf")
names(world) <- tolower(names(world))
world <- world[world$adm0_a3 != "ATA", ]                 # drop Antarctica

status_by_code <- unique(merge(lookup, ec[, c("Region", "Current_class_proposed")],
                               by = "Region")[, c("adm0_a3", "Current_class_proposed")])
status_by_code <- status_by_code[status_by_code$adm0_a3 %in% lookup$adm0_a3[lookup$Region %in% ipd$Region], ]
# One status per map polygon (France and French Guiana are both non-endemic)
if (any(duplicated(status_by_code$adm0_a3))) {
  stop("Conflicting status for one map polygon: ",
       paste(status_by_code$adm0_a3[duplicated(status_by_code$adm0_a3)], collapse = ", "))
}
# Somaliland is a separate polygon in Natural Earth; it is part of Somalia here
if ("SOM" %in% status_by_code$adm0_a3) {
  status_by_code <- rbind(status_by_code,
                          data.frame(adm0_a3 = "SOL",
                                     Current_class_proposed =
                                       status_by_code$Current_class_proposed[status_by_code$adm0_a3 == "SOM"]))
}
world <- merge(world, status_by_code, by = "adm0_a3", all.x = TRUE)
world$fill_class <- ifelse(is.na(world$Current_class_proposed), "No included patients",
                           world$Current_class_proposed)
world$fill_class <- factor(world$fill_class,
                           levels = c("Endemic", "Non-endemic", "No included patients"))
not_found <- setdiff(status_by_code$adm0_a3, world$adm0_a3)
if (length(not_found)) warning("Map codes not found in Natural Earth: ",
                               paste(not_found, collapse = ", "))

bubbles <- merge(counts, lookup, by = "Region")
bubbles <- st_as_sf(bubbles, coords = c("lon", "lat"), crs = 4326)
bubbles <- bubbles[order(-bubbles$n), ]                   # small bubbles drawn on top

# Zoom area for panel B: Mediterranean and Middle East
ZOOM <- c(xmin = -12, xmax = 62, ymin = 14, ymax = 50)
box_x <- c(seq(ZOOM["xmin"], ZOOM["xmax"], length.out = 50), rep(ZOOM["xmax"], 50),
           seq(ZOOM["xmax"], ZOOM["xmin"], length.out = 50), rep(ZOOM["xmin"], 50))
box_y <- c(rep(ZOOM["ymin"], 50), seq(ZOOM["ymin"], ZOOM["ymax"], length.out = 50),
           rep(ZOOM["ymax"], 50), seq(ZOOM["ymax"], ZOOM["ymin"], length.out = 50))
zoom_box <- st_sfc(st_polygon(list(cbind(c(box_x, box_x[1]), c(box_y, box_y[1])))),
                   crs = 4326)   # many vertices so the box bends correctly in Robinson
in_zoom <- sum(bubbles$n[st_coordinates(bubbles)[, 1] >= ZOOM["xmin"] &
                         st_coordinates(bubbles)[, 1] <= ZOOM["xmax"] &
                         st_coordinates(bubbles)[, 2] >= ZOOM["ymin"] &
                         st_coordinates(bubbles)[, 2] <= ZOOM["ymax"]])
cat("Patients inside the zoom area (panel B):", in_zoom, "of", sum(bubbles$n), "\n")

# -----------------------------------------------------------------------------
# 4. Plot: A = world (Robinson), B = Mediterranean and Middle East
# -----------------------------------------------------------------------------
library(patchwork)
FILL <- c("Endemic" = "#C9733A", "Non-endemic" = "#2F7DBF",
          "No included patients" = "#E6E6E6")
size_breaks <- c(1, 5, 20, 60)
size_scale  <- scale_size_area(max_size = 13, breaks = size_breaks,
                               limits = c(0, max(bubbles$n)), name = "Patients")

base_layers <- list(
  geom_sf(data = world, aes(fill = fill_class), colour = "white", linewidth = 0.15),
  geom_sf(data = bubbles, aes(size = n), shape = 21, fill = "#1F1F1F",
          colour = "white", stroke = 0.4, alpha = 0.7),
  scale_fill_manual(values = FILL, name = "Current status", drop = FALSE),
  size_scale,
  guides(fill = guide_legend(order = 1, override.aes = list(size = 0.5)),
         size = guide_legend(order = 2, override.aes = list(fill = "#1F1F1F"))),
  theme_void(base_size = 9),
  theme(legend.title = element_text(face = "bold"),
        plot.margin = margin(4, 4, 4, 4))
)

p_world <- ggplot() + base_layers +
  geom_sf(data = zoom_box, fill = NA, colour = "grey25", linewidth = 0.35, linetype = "22") +
  coord_sf(crs = "+proj=robin", expand = FALSE)

# Country labels for panel B: hand-placed so they do not overlap; a thin leader
# line joins small countries to their label
zoom_labels <- read.csv(text = "
Region,lab_lon,lab_lat,hjust,leader
Spain,-3.7,37.4,0.5,N
Portugal,-8.4,43.4,0.5,N
France,2.4,48.2,0.5,N
Czech Republic,18.6,48.2,0,Y
Italy,9.6,41.0,1,Y
Bosnia and Herzegovina,19.8,46.0,0,Y
Greece,21.0,36.2,0.5,N
Turkey,35.0,42.6,0.5,N
Syria,41.2,36.6,0,Y
Lebanon,31.6,34.8,1,Y
Israel,31.6,31.8,1,Y
Jordan,39.0,29.6,0,Y
Kuwait,44.6,31.4,1,Y
Qatar,53.2,27.4,0,Y
Oman,55.0,17.6,1,Y
Saudi Arabia,45.0,21.6,0.5,N
Iran,53.0,30.2,0.5,N
Tunisia,12.2,33.0,0,Y
Algeria,2.6,26.4,0.5,N
Morocco,-6.5,30.2,0.5,N
Armenia,49.5,39.0,0,Y
Russia,49.5,45.5,0,Y
", stringsAsFactors = FALSE, strip.white = TRUE)
zoom_labels <- merge(zoom_labels, merge(counts, lookup, by = "Region"), by = "Region")
zoom_labels <- zoom_labels[zoom_labels$lon >= ZOOM["xmin"] & zoom_labels$lon <= ZOOM["xmax"] &
                           zoom_labels$lat >= ZOOM["ymin"] & zoom_labels$lat <= ZOOM["ymax"], ]
zoom_labels$label <- paste0(ifelse(zoom_labels$Region == "Bosnia and Herzegovina",
                                   "Bosnia and Herz.", zoom_labels$Region),
                            " (", zoom_labels$n, ")")
unlabelled <- setdiff(
  bubbles$Region[st_coordinates(bubbles)[, 1] >= ZOOM["xmin"] & st_coordinates(bubbles)[, 1] <= ZOOM["xmax"] &
                 st_coordinates(bubbles)[, 2] >= ZOOM["ymin"] & st_coordinates(bubbles)[, 2] <= ZOOM["ymax"]],
  zoom_labels$Region)
if (length(unlabelled)) message("No label position set for: ", paste(unlabelled, collapse = ", "))
leaders <- zoom_labels[zoom_labels$leader == "Y", ]

p_zoom <- ggplot() + base_layers +
  geom_segment(data = leaders, aes(x = lon, y = lat, xend = lab_lon, yend = lab_lat),
               colour = "grey25", linewidth = 0.25) +
  geom_text(data = zoom_labels, aes(x = lab_lon, y = lab_lat, label = label, hjust = hjust),
            size = 2.5, colour = "grey10") +
  coord_sf(xlim = ZOOM[c("xmin", "xmax")], ylim = ZOOM[c("ymin", "ymax")],
           crs = st_crs(4326), expand = FALSE) +
  theme(panel.border = element_rect(colour = "grey25", fill = NA, linewidth = 0.35,
                                    linetype = "22"))

p <- (p_world / p_zoom) +
  plot_layout(heights = c(1.15, 1), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom", legend.box = "horizontal",
        plot.tag = element_text(face = "bold", size = 11))

print(p)

out_dir <- file.path("outputs", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
w <- 10; h <- 10
ggsave(file.path(out_dir, "Fig2_world_map.pdf"),  p, width = w, height = h)
ggsave(file.path(out_dir, "Fig2_world_map.png"),  p, width = w, height = h, dpi = 600)
ggsave(file.path(out_dir, "Fig2_world_map.tiff"), p, width = w, height = h, dpi = 600,
       compression = "lzw")
cat("Saved to", normalizePath(out_dir), "\n")
