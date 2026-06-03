# ============================================================
# PPI Network Analysis — BaP → Lung Cancer (30 intersection genes)
# STRING v12, Cytoscape-compatible export
# ============================================================

library(httr)
library(jsonlite)
library(igraph)
library(ggraph)
library(tidyverse)

setwd("C:/Users/1/Desktop/17(1)")

# ---- 1. Load intersection genes ----
genes <- read.table("interGenes.txt", header = FALSE, stringsAsFactors = FALSE)[,1]
genes <- genes[genes != ""]  # remove empty
cat(sprintf("Input genes: %d\n", length(genes)))
cat(paste(genes, collapse = ", "), "\n")

# ---- 2. STRING API call ----
string_api <- function(genes, species = 9606, score_threshold = 400) {
  url <- "https://string-db.org/api/json/network"
  params <- list(
    identifiers = paste(genes, collapse = "%0d"),
    species = species,
    required_score = score_threshold,
    network_type = "functional"
  )
  resp <- GET(url, query = params)
  stop_for_status(resp)
  fromJSON(content(resp, as = "text", encoding = "UTF-8"))
}

cat("\nFetching PPI from STRING v12...\n")
ppi <- string_api(genes, score_threshold = 400)
cat(sprintf("Retrieved %d interactions among %d proteins\n", nrow(ppi),
            length(unique(c(ppi$preferredName_A, ppi$preferredName_B)))))

# ---- 3. Build igraph ----
edges <- ppi[, c("preferredName_A", "preferredName_B", "score")]
colnames(edges) <- c("from", "to", "combined_score")
edges$combined_score <- as.numeric(edges$combined_score) / 1000  # scale to 0-1

g <- graph_from_data_frame(edges, directed = FALSE)
E(g)$weight <- edges$combined_score
g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

cat(sprintf("Final network: %d nodes, %d edges\n", vcount(g), ecount(g)))

# ---- 4. Topological analysis ----
V(g)$degree        <- degree(g)
V(g)$betweenness   <- betweenness(g, normalized = TRUE)
V(g)$closeness     <- closeness(g, normalized = TRUE)
V(g)$eigenvector   <- eigen_centrality(g)$vector
V(g)$clustering    <- transitivity(g, type = "local")

# Identify hub genes (top 10 by degree)
node_stats <- data.frame(
  Gene = V(g)$name,
  Degree = V(g)$degree,
  Betweenness = round(V(g)$betweenness, 4),
  Closeness = round(V(g)$closeness, 4),
  Eigenvector = round(V(g)$eigenvector, 4),
  ClusteringCoeff = round(V(g)$clustering, 4),
  stringsAsFactors = FALSE
)
node_stats <- node_stats[order(-node_stats$Degree), ]
rownames(node_stats) <- NULL

cat("\n========== Top 10 Hub Genes (by Degree) ==========\n")
print(head(node_stats, 10))

# ---- 5. Export node/edge tables (Cytoscape compatible) ----
edge_out <- data.frame(
  source = edges$from,
  target = edges$to,
  combined_score = edges$combined_score,
  stringsAsFactors = FALSE
)

write.table(node_stats, "PPI_nodes.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(edge_out, "PPI_edges.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(node_stats, "PPI_nodes.csv", row.names = FALSE)
write.csv(edge_out, "PPI_edges.csv", row.names = FALSE)
cat("\nNode/edge tables exported: PPI_nodes.txt, PPI_edges.txt\n")

# ---- 6. Visualization ----
# Color gradient by degree
V(g)$color <- colorRampPalette(c("#f0f0f0", "#2171b5", "#08306b"))(max(V(g)$degree) + 1)[V(g)$degree + 1]
V(g)$size  <- scales::rescale(V(g)$degree, to = c(3, 12))
V(g)$label <- ifelse(V(g)$degree >= quantile(V(g)$degree, 0.75), V(g)$name, "")

set.seed(42)
pdf("PPI_network.pdf", width = 10, height = 9)

ggraph(g, layout = "fr") +
  geom_edge_link(aes(alpha = weight), color = "#bdbdbd", show.legend = FALSE) +
  geom_node_point(aes(size = degree, fill = degree), shape = 21, color = "white", stroke = 0.3) +
  geom_node_text(aes(label = label), size = 3.2, repel = TRUE,
                 max.overlaps = 30, box.padding = 0.6) +
  scale_fill_gradient(low = "#deebf7", high = "#08306b", name = "Degree") +
  scale_size_continuous(range = c(3, 12), name = "Degree") +
  labs(
    title = "PPI Network: BaP Target Genes in Lung Cancer",
    subtitle = paste0(vcount(g), " nodes, ", ecount(g), " edges | STRING v12 (score ≥ 0.4)"),
    caption = "Node size/color = degree; Labeled = top 25% hub genes"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
    plot.caption = element_text(hjust = 0.5, size = 8, color = "grey60"),
    legend.position = "right"
  )

dev.off()
cat("PPI network PDF saved: PPI_network.pdf\n")

# ---- 7. Hub gene barplot ----
top15 <- head(node_stats, 15)
pdf("PPI_hub_barplot.pdf", width = 10, height = 6)

ggplot(top15, aes(x = reorder(Gene, Degree), y = Degree, fill = Degree)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = Degree), hjust = -0.3, size = 3.5) +
  scale_fill_gradient(low = "#6baed6", high = "#08306b") +
  coord_flip() +
  labs(
    title = "Top 15 Hub Genes in PPI Network",
    subtitle = "BaP → Lung Cancer Intersection Genes",
    x = "", y = "Degree Centrality"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

dev.off()
cat("Hub gene barplot saved: PPI_hub_barplot.pdf\n")

# ---- 8. Summary ----
cat("\n========== PPI Analysis Complete ==========\n")
cat(sprintf("Input genes: %d\n", length(genes)))
cat(sprintf("PPI nodes: %d, edges: %d\n", vcount(g), ecount(g)))
cat(sprintf("Network density: %.4f\n", edge_density(g)))
cat(sprintf("Mean degree: %.1f\n", mean(V(g)$degree)))
cat(sprintf("Top hub: %s (degree=%d)\n", node_stats$Gene[1], node_stats$Degree[1]))
