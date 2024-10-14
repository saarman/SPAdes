# Load required libraries
library(ggplot2)
library(data.table)

# Define the input file path
input_file <- "/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/output/all_DB_lin_clu.tsv"

# Define the output directory
output_dir <- "/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/output/"  # Update this path to your desired output directory

# Read the TSV file
data <- fread(input_file, header = FALSE)

# Assign column names (assuming the format is representative sequence and cluster members)
colnames(data) <- c("Representative", "ClusterMember")

# Create a data frame for plotting
plot_data <- data.table::melt(data, id.vars = "Representative", variable.name = "Type", value.name = "Sequence")

# Generate the plot
p <- ggplot(plot_data, aes(x = Representative, y = Sequence)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Clustering Results", x = "Representative Sequence", y = "Cluster Members")

# Save the plot to a file in the specified output directory
output_file <- file.path(output_dir, "clustering_plot.png")
ggsave(output_file, plot = p, width = 10, height = 8)
