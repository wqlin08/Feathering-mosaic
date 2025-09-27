# Feathering-mosaic
A feathering method for tif data mosaicking to eliminate splicing artifacts
# Grid-Based Raster Mosaicking with Feathering

This Python script provides a high-performance solution for mosaicking a large number of raster tiles (e.g., GeoTIFFs from drone surveys) into a seamless, grid-aligned dataset. It is designed to handle massive datasets containing tens or hundreds of thousands of files efficiently.

## Key Features

-   **Efficient Spatial Indexing**: Creates a spatial index of all input tiles for rapid querying, drastically speeding up processing for large areas.
-   **Seamless Feathering**: Implements a weighted average algorithm to blend overlapping tiles, eliminating sharp seamlines and creating a visually continuous mosaic.
-   **Grid-Based Processing**: Mosaics tiles based on a user-provided grid shapefile, ensuring standardized outputs.
-   **Configurable Overlap**: Adds a specified overlap (buffer) to each output grid tile, which is useful for subsequent analysis or further mosaicking.
-   **Automated Workflow**: Scans multiple input directories and processes all grids automatically. Skips existing files for easy resumption of interrupted jobs.

## Core Algorithms

### 1. Spatial Indexing (R-Tree)

To avoid iterating through every single input file for each grid cell, the script first builds a spatial index.

-   **Process**: It reads the geographic bounds of every source `.tif` file and stores them in a GeoPackage file using GeoPandas.
-   **Benefit**: For large datasets, this is the most critical performance optimization.

### 2. Seamless Mosaicking (Feathering)

To eliminate seams where tiles overlap, a feathering technique based on a weighted average is used.

-   **Process**:
    1.  A "weight map" is generated for each tile, where pixels at the center have a weight of 1.0, and pixels near the edges have weights that fade smoothly to 0.
    2.  When mosaicking, two buffers are accumulated: one for the sum of `pixel_value * weight` and another for the sum of `weight`.
    3.  The final pixel value is calculated as `total_weighted_value / total_weight`.
-   **Benefit**: This ensures a smooth, gradual transition in overlapping areas instead of an abrupt line, resulting in a high-quality, seamless final product.

## 3. Customization and Scalability

This script is optimized for large-scale, processing where the source dataset is too large to be handled as a single file. The grid-based approach is essential when dealing with massive datasets (e.g., hundreds of thousands of tiles) as it breaks the problem down into manageable chunks.

However, it can be easily adapted for simpler tasks.

### For a Single, Large Mosaic (No Grids)

If you have a smaller number of tiles and want to create one single output GeoTIFF instead of a grid of tiles, you can modify the main execution block. You can remove the grid logic and the related parts of the multipath input and only keep the related parts of the feathered mosaic.
