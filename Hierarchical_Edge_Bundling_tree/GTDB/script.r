setwd('.')

# Switch to R
library(tidyverse)
library(viridis)
library(patchwork)
library(ggraph)
library(igraph)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(rlang)


##### Defining variables #####

# Define total random postive and negative estimations to be displayed in tree
n=50 
# Define the level of taxonomy (class, genus, family, etc))
taxonomy_level_column<-"family" # Can change to any other taxonomy levels such as family
# build the dynamic column names
taxon_from_col <- sym(paste0(taxonomy_level_column, "_from"))
taxon_to_col   <- sym(paste0(taxonomy_level_column, "_to"))
# Define file path to genomic dnds values
file_path <- "/data/jzr5814/sourmash_dnds_estimation/for_jinglin/fmh_omega_7.csv"

###############################
###############################


##### Step 1: Read in the raw taxonomy files for both Archaea and Bacteria #####
# Taxonomy for Archaea (Read in the raw taxonomy file)
archaea_taxonomy_raw <- read.csv('GTDB_taxonomy/ar53_taxonomy_r214.tsv', sep = '\t', header = FALSE)
colnames(archaea_taxonomy_raw) <- c('Assembly ID', 'lineage')

# Taxonomy for Bacteria (Read in the raw taxonomy file)
bacteria_taxonomy_raw <- read.csv('GTDB_taxonomy/bac120_taxonomy_r214.tsv', sep = '\t', header = FALSE)
colnames(bacteria_taxonomy_raw) <- c('Assembly ID', 'lineage')

###############################
###############################


##### Step 2: Clean the archaea taxonomy raw data #####
archaea_taxonomy <- archaea_taxonomy_raw %>%
  mutate(`Assembly ID` = gsub("^.{2}_", "", `Assembly ID`)) %>%  # Remove first two characters and underscore
  mutate(root = "Root") %>%  # Add the "root" column with all values as "Root"
  relocate(root, .after = `Assembly ID`) %>%  # Place the "root" column after "Assembly ID"
  separate(lineage, into = c("superkingdom", "phylum", "class", "order", "family", "genus", "species"),
           sep = ";", remove = TRUE) %>%  # Split lineage into taxonomy columns
  mutate(across(superkingdom:species, ~ gsub("^[dpcfogs]__", "", .)))


# Clean the bacteria taxonomy raw data
bacteria_taxonomy <- bacteria_taxonomy_raw %>%
  mutate(`Assembly ID` = gsub("^.{2}_", "", `Assembly ID`)) %>%  # Remove first two characters and underscore
  mutate(root = "Root") %>%  # Add the "root" column with all values as "Root"
  relocate(root, .after = `Assembly ID`) %>%  # Place the "root" column after "Assembly ID"
  separate(lineage, into = c("superkingdom", "phylum", "class", "order", "family", "genus", "species"),
           sep = ";", remove = TRUE) %>%  # Split lineage into taxonomy columns
  mutate(across(superkingdom:species, ~ gsub("^[dpcfogs]__", "", .)))

# Append the archaea and bacteria taxonomy data
lineage_information_new <- rbind(archaea_taxonomy, bacteria_taxonomy)

#colnames(lineage_information_new)
#typeof(lineage_information_new)

# Export to CSV
#write.csv(lineage_information_new, file = "lineage_information.csv", row.names = FALSE, quote = FALSE)

###############################
###############################


##### Step 3: Clean the 'query_name_x' and 'match_name_x' columns in the 'dnds_results' dataset #####

# Read the CSV file
dnds_results <- read.csv(file_path, sep = ",", header = TRUE)

# Modify 'query_name_x' and 'match_name_x' columns to remove the first two characters and underscore
dnds_results$query_name <- sapply(strsplit(as.character(dnds_results$query_name), "_"), function(x) paste(x[-1], collapse = "_"))
dnds_results$match_name <- sapply(strsplit(as.character(dnds_results$match_name), "_"), function(x) paste(x[-1], collapse = "_"))

# Filter rows where 'query_name_x' and 'match_name_x' are in the column "Assembly ID" of lineage_information_new
dnds_results_modified <- dnds_results[
  dnds_results$query_name %in% lineage_information_new$`Assembly ID` &
  dnds_results$match_name %in% lineage_information_new$`Assembly ID`,
]

#colnames(lineage_information_new)
#typeof(lineage_information_new)

# Select specific columns: 'query_name_x', 'match_name_x', and 'dN/dS'
dnds_results_modified <- dnds_results_modified[, c("query_name", "match_name", "dN.dS")]

# Export to CSV
#write.csv(dnds_results_modified, file = "dnds_results_modified.csv", row.names = FALSE, quote = FALSE)

###############################
###############################


##### Step 4: Clean taxonomy and dnds datasets, selecting columns that are needed #####
#taxonomy <- read.csv('lineage_information_jzr.csv')
#colnames(taxonomy)
#typeof(taxonomy)
#dnds <- read.csv('fmh_omega_7_short_jzr.csv') #this can be erased # deprecated

dnds <- dnds_results_modified %>%                    # change dnds with dnds_results_modified
  filter(!is.na(dN.dS)) %>%
  rename(from = query_name, to = match_name, dndsvalue = dN.dS) %>%
  select(from, to, dndsvalue)


# Keep only Assembly.ID present in dnds$query_name_x or dnds$match_name_x
taxonomy_selected <- lineage_information_new %>% # change taxononmy with lineage_information_new
  filter(`Assembly ID` %in% unique(c(dnds$from, dnds$to)))

###############################
###############################


##### Step 5: Check for mutual pairs such as from genome A to genome B and from genome B to genome A in dnds dataset #####
mutual_check <- merge(dnds, dnds, by.x = c("from", "to"), by.y = c("to", "from"))

# Output the result
if (nrow(mutual_check) > 0) {
  print("Mutual connections found:")
  print(mutual_check)
} else {
  print("No mutual connections found.")
}

###############################
###############################


##### Step 6: # Find out all genomes ids that are in the isolated pair #####

# The isolated pair here means each genome in this pair only connects to other genomes that only exists in this pair

# Create an undirected graph from the edge list
g <- graph_from_data_frame(dnds, directed = FALSE)

# Find connected components
comps <- components(g)

# Get component sizes
component_sizes <- comps$csize

# Get the membership vector
membership <- comps$membership

# Find components of size 2 (which connects with one and only one other genome)
components_of_size_two <- which(component_sizes == 2)

# Initialize vector to store genome IDs in isolated pairs
isolated_genomes <- c()

# Loop through each component of size 2
for (comp in components_of_size_two) {
  # Get the nodes in this component
  nodes_in_comp <- names(membership[membership == comp])
  # Add them to isolated_genomes
  isolated_genomes <- c(isolated_genomes, nodes_in_comp)
}

# Remove duplicate entries
isolated_genomes <- unique(isolated_genomes)

###############################
###############################


##### deprecated #####
# # Step 7: # Based on the dnds dataset, removing all species that contained isolated_genomes ID in Step 6 
# since we do not want to show gemones that connects only one another genome in the visualization
# filter dnds to keep the visulation clean
# connect <- dnds[!(dnds$from %in% isolated_genomes), ] %>%
  # filter((dndsvalue >= 0.7 & dndsvalue <= 0.9 | dndsvalue >= 1.1 & dndsvalue <= 3))
  # filter(dndsvalue > 0.8 & dndsvalue < 1)
  # filter(dndsvalue > 1)
  # filter(dndsvalue >= 1.1 & dndsvalue <= 2)

###############################
###############################


##### Step 7: Let's randomly sort and randomly choose so that the dataset isnt so big #####

# Filter to exclude isolated genomes
connect_temp <- dnds[!(dnds$from %in% isolated_genomes), ]

# Randomly select 50 rows for each dnds range
n=50
set1 <- connect_temp %>%
  filter(dndsvalue >= 0.4, dndsvalue <= 0.9) %>%
  sample_n(n)

set2 <- connect_temp %>%
  filter(dndsvalue >= 1.1, dndsvalue <= 2) %>%
  sample_n(n)

# Combine the two sets
connect <- bind_rows(set1, set2)
connect_csv_file <- paste0("connect_", taxonomy_level_column, ".csv")
readr::write_csv(connect, connect_csv_file)

###############################
###############################


##### Step 8: Merge taxonomy information with the 'connect' dataset to #####

# ensure both genome A and genome B have taxonomy information for a given taxonomic level such as genus

tax_cols <- c("root","superkingdom","phylum","class","order","family","genus","species")  

tax_to <- taxonomy_selected %>%
  rename(Assembly.ID_to = `Assembly ID`) %>%                                         
  rename_with(~ paste0(.x, "_to"), all_of(tax_cols))                                

tax_from <- taxonomy_selected %>%
  rename(Assembly.ID_from = `Assembly ID`) %>%                                       
  rename_with(~ paste0(.x, "_from"), all_of(tax_cols))                              

# join taxonomy for BOTH ends so *_from columns exist
connect_with_taxonomy <- connect %>%
  left_join(tax_from, by = c("from" = "Assembly.ID_from")) %>%                      
  left_join(tax_to,   by = c("to"   = "Assembly.ID_to")) %>%                        
  mutate(
    taxon_from = !!taxon_from_col,                                                  
    taxon_to   = !!taxon_to_col                                                     
  )

# Now your filtering works because superkingdom_from/_to and taxon_from/_to exist
connect_with_taxonomy_updated <- connect_with_taxonomy %>%
  filter(
    taxon_from != "" & !is.na(taxon_from),
    taxon_to   != "" & !is.na(taxon_to),
    superkingdom_from == superkingdom_to
  )

cat("Number of valid edges (with positions):", nrow(connect_with_taxonomy_updated), "\n")

connect_with_taxonomy_updated_same_genus <- connect_with_taxonomy_updated %>%
  filter(taxon_from == taxon_to) %>%
  arrange(taxon_from)


# write to csv
# write.csv(connect_with_taxonomy_updated_same_genus, file = "dnds_0.8_1_same_genus_jzr.csv", row.names = FALSE, quote = FALSE)

###############################
###############################


##### Step 9: Clean the taxonomy_selected df to keep only the nodes (valid nodes) that are present #####
# in the filtered connect df in Step 8
taxonomy_selected_final <- taxonomy_selected[taxonomy_selected$`Assembly ID` %in% unique(c(connect_with_taxonomy_updated$from, connect_with_taxonomy_updated$to)), ]
cat("Number of valid nodes:", nrow(taxonomy_selected_final), "\n")

###############################
###############################


##### Step 10: Create hierarcy and vertices for hierarchical edge bundling visualization #####

# Add a root column to the taxonomy data, reshape the data into long format
taxonomy_long <- taxonomy_selected_final %>%
  select(`Assembly ID`, root, superkingdom, phylum, class, order, family, genus, species) %>%
  mutate(root = if_else(root == "cellular organisms", "Root", root)) %>%
      # Create hierarchical relationships
  pivot_longer(cols = c(root, superkingdom, phylum, class, order, family, genus, species, `Assembly ID`),
               names_to = "level",
               values_to = "name") %>%
  mutate(genome_id = rep(taxonomy_selected_final$`Assembly ID`, each = 9)) # 9 is the numbers of columns in taxonomy_selected

  # Create hierarchical relationships
hierarchy <- taxonomy_long %>%

  # sorting taxomomic level within each genome ("genome_id")
  arrange(genome_id, match(level, c("root", "superkingdom", "phylum", "class", "order", "family", "genus", "species", "Assembly ID"))) %>%
  group_by(genome_id) %>%

  # from = lag(name): Take the value from the previous row in the column "name" and assigns it to a new column called from, 
  # to = name: Set the value in column "to" the value in the name column.
  mutate(from = lag(name), to = name)  %>%
  filter(!is.na(from)) %>%
  select(from, to, genome_id) %>%
  distinct() %>%
  ungroup() %>%
  distinct(from, to)

vertices <- data.frame(name = unique(c(as.character(hierarchy$from), as.character(hierarchy$to))))


# Create the graph object
vertices <- data.frame(name = unique(c(as.character(hierarchy$from), as.character(hierarchy$to))), value = runif(nrow(vertices))) 
mygraph <- graph_from_data_frame(hierarchy, vertices = vertices)


from <- match(connect_with_taxonomy_updated$from, vertices$name)
to <- match(connect_with_taxonomy_updated$to, vertices$name)

###############################
###############################


##### Step 11: Retain layout data for only the nodes that are connected by valid edges #####
# create the layout so we can feed to to get_con and manipulate the data frame directly for plotting
layout_data <- create_layout(mygraph, layout = 'dendrogram', circular = TRUE)

# # 'from' and 'to' are vectors of node indices or names
conn_function <- get_con(from = from, to = to)
conn_data_frame <- conn_function(layout_data)
edges <- conn_data_frame


# # # Add connection values to `conn_data_frame`. This is one of the key steps of getting the connections/dNdS values to match up with the actual plot
conn_data_frame$con.value <- connect_with_taxonomy_updated$dndsvalue[conn_data_frame$con.id]

# Remove invalid connections where genomes are not in connect_new
conn_data_frame <- conn_data_frame[!is.na(conn_data_frame$con.value), ]

sum(conn_data_frame$leaf == 'TRUE')

###############################
###############################


##### Step 12: Determine the coordinates of the genus label #####
genome_coordinate<-unique(conn_data_frame[conn_data_frame$leaf == TRUE, c('name','x','y')])

genome_coordinate_info <- merge(
  genome_coordinate, 
  taxonomy_selected_final[, c("Assembly ID", "phylum", "class", "order", "family", "genus", "species")],
  by.x = "name", 
  by.y = "Assembly ID", 
  all.x = TRUE
)

# Function to calculate center coordinates for each genus
calculate_taxon_centers <- function(data, taxonomy_level_column, x_column, y_column) {
  # Group by genus and calculate the mean x and y coordinates
  centers <- data %>%
    group_by(!!sym(taxonomy_level_column)) %>%
    summarise(
      center_x = mean(!!sym(x_column)),
      center_y = mean(!!sym(y_column)),
      .groups = 'drop'
    )
  return(centers)
}

# Dynamically create the name for taxon_centers
centers_name <- paste0("taxon", "_centers")

# Calculate centers and assign to dynamically named variable
assign(centers_name, calculate_taxon_centers(genome_coordinate_info, 
                                             taxonomy_level_column = taxonomy_level_column, 
                                             x_column = "x", 
                                             y_column = "y"))


cat("Number of unique", taxonomy_level_column, "groups:", nrow(taxon_centers), "\n")

# node labels
taxon_centers <- taxon_centers %>%
  mutate(
    # raw angle in degrees from center (0 deg = +x axis, increases ccw)
    angle = atan2(center_y, center_x) * 180 / pi,
    angle = if_else(angle < 0, angle + 360, angle),   

    # flip labels on the left side so they are upright
    hjust = if_else(angle > 90 & angle < 270, 1, 0),
    angle_text = if_else(hjust == 1, angle + 180, angle),

    # pre-scale label radius so you can tweak outward offset in one place
    label_x = center_x * 1.15,
    label_y = center_y * 1.15
  )

###############################
###############################


##### Step 13: Generate a graph based on genomes and their dN/dS values, labelling each dot on the circle represents a genome #####

p <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_conn_bundle(
    data = conn_data_frame, 
    alpha = 0.8, 
    width = 0.9, 
    aes(edge_colour = con.value)
  ) +
geom_text(
    data = taxon_centers,
    aes(
      x = label_x,                
      y = label_y,                 
      label = !!sym(taxonomy_level_column),
      angle = angle_text,          
      hjust = hjust         
    ),
    size = 2,
    alpha = 1,
    check_overlap = TRUE         
  ) +
  geom_node_point(
    aes(
      filter = leaf,
      x = x * 1.07,
      y = y * 1.07
    ),
    alpha = 0.2
  ) +
  scale_edge_colour_gradient2(
    low = "blue",
    mid = "white",
    high = "darkorange",
    midpoint = 1,
    name = "FracMinHash dN/dS"
  ) +
  scale_size_continuous(range = c(1, 3)) +
  theme_void() +
  theme(
    plot.margin = margin(5, 12, 5, 12),   
    legend.position = "bottom",             
    legend.justification = "right",         
    legend.box.just = "right",
    legend.direction = "horizontal",
    panel.background = element_rect(fill = "white", colour = NA),  
    plot.background  = element_rect(fill = "white", colour = NA)    
  ) +
  expand_limits(x = c(-1.8, 1.8), y = c(-1.8, 1.8))

###############################
###############################


##### Step 14: Save as a png #####
png_file <- paste0("output_", taxonomy_level_column, ".png")
ggsave(png_file, plot = p, width = 10, height = 10, dpi = 300)

