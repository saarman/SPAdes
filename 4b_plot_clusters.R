# Load required libraries
library(ggplot2)
library(data.table)
library(igraph)
library(ggraph)

# Define the input file path
input_file <- "/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/output/all_DB_clu.tsv"

# Define the output directory
output_dir <- "/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/output/"  # Update this path to your desired output directory

# Read the TSV file
data <- fread(input_file, header = FALSE)

# Assign column names (assuming the format is representative sequence and cluster members)
colnames(data) <- c("Representative", "ClusterMember")

# Create an edge list for the network graph
edge_list <- data[, .(Representative, ClusterMember)]

# Create an igraph object from the edge list
network_graph <- graph_from_data_frame(edge_list, directed = FALSE)

# Generate layout coordinates
layout <- create_layout(network_graph, layout = "fr")

# Plot the network without any labeling
p <- ggraph(layout) +  # Use the layout with coordinates
  geom_edge_link(col = "black", width = 0.5) +  # Edges in black
  geom_node_point(size = 2) +  # Size of nodes
  theme_minimal() +  # Minimal theme
  ggtitle("Clustering Network")

# Save the plot to a file in the specified output directory
output_file <- file.path(output_dir, "network_plot.png")
ggsave(output_file, plot = p, width = 10, height = 7)
