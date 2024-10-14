# Load required libraries
library(ggplot2)
library(data.table)
library(igraph)
library(ggraph)

# Define the input file path
input_file <- "/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/output/all_DB_lin_clu.tsv"

# Define the output directory
output_dir <- "/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/output/"  # Update this path to your desired output directory

# Read the TSV file
data <- fread(input_file, header = FALSE)

# Assign column names (assuming the format is representative sequence and cluster members)
colnames(data) <- c("Representative", "ClusterMember")

# Convert Representative to a factor
data$Representative <- as.factor(data$Representative)

# Create an edge list for the network graph
edge_list <- data[, .(Representative, ClusterMember)]

# Create an igraph object from the edge list
network_graph <- graph_from_data_frame(edge_list, directed = FALSE)

# Print a summary of the graph to verify it loaded correctly
print(summary(network_graph))
print(network_graph, e = TRUE, v = TRUE)

# Prepare the edge data for plotting
edge_data <- as_data_frame(network_graph, what = "edges")
edge_data$Representative <- edge_list$Representative

# Prepare the vertex data for plotting
vertex_data <- as_data_frame(network_graph, what = "vertices")

# Plot the network with edges colored by representative sequence
p <- ggraph(network_graph, layout = "fr") +  # "fr" is the Fruchterman-Reingold layout
  geom_edge_link(aes(color = Representative), data = edge_data, width = 1) +  # Color edges by representative sequence
  geom_node_point(aes(), data = vertex_data, size = 2) +  # Size of nodes
  geom_node_text(aes(label = name), data = vertex_data, vjust = 1, hjust = 1) +  # Label nodes
  scale_edge_color_manual(values = rainbow(length(unique(data$Representative)))) +  # Color palette
  theme_minimal() +  # Minimal theme
  ggtitle("Clustering Network")

# Save the plot to a file in the specified output directory
output_file <- file.path(output_dir, "network_plot.png")
ggsave(output_file, plot = p, width = 10, height = 7)
