FROM rocker/r-ver:4.6.0

# System libraries that many CRAN/Bioc packages need to compile.
# Add/remove based on your actual packages.
RUN apt-get update && apt-get install -y --no-install-recommends \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libfontconfig1-dev \
    && rm -rf /var/lib/apt/lists/*

# --- CRAN packages ---
# Rocker images use a pinned RSPM/Posit snapshot, so versions are reproducible.
RUN R -q -e 'install.packages(c("dplyr", "tidyverse", "ggplot2", "ggrepel", "DT", "openxlsx", "rmarkdown"))'

# --- Bioconductor packages ---
RUN R -q -e 'install.packages("BiocManager")' \
 && R -q -e 'BiocManager::install(c("PCAtools", "Biostrings"), update = FALSE, ask = FALSE)'

# --- Local package from a zip archive ---
# Copy the archive into the image, then install it.
COPY MSstatsAlt_4.17.6.tar.gz /tmp/MSstatsAlt_4.17.6.tar.gz
#RUN R -q -e 'install.packages("/tmp/MSstatsAlt_4.17.6.tar.gz", repos = NULL, type = "source")' 
RUN R -q -e ' \
    install.packages("remotes"); \
    remotes::install_local("/tmp/MSstatsAlt_4.17.6.tar.gz", dependencies = TRUE); \
    if (!requireNamespace("MSstatsAlt", quietly = TRUE)) stop("install failed") \
    '