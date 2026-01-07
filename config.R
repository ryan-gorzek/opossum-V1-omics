#### config.R
## 1. Read YAML config
## 2. Install/load comparatome
## 3. Build dir.list with absolute paths
## 4. Create directories if needed
## 5. (Optional) Download processed GEO data

if (!requireNamespace("yaml", quietly = TRUE)) {
  install.packages("yaml")
}
library(yaml)

# Read YAML
cfg <- yaml::read_yaml("config.yml")   # assumes config.R is run from project root or use here::here()

# Basic scalar config
tools.path <- cfg$tools_path
root.path  <- cfg$root_path

# Convert relative paths under cfg$dir into absolute paths under root.path
make_abs_paths <- function(x, root) {
  add_trailing_slash <- function(p) {
    if (is.character(p) && length(p) == 1L && nzchar(p) && !grepl("/$", p)) {
      paste0(p, "/")
    } else {
      p
    }
  }
  
  if (is.list(x)) {
    lapply(x, make_abs_paths, root = root)
  } else if (is.character(x) && length(x) == 1L) {
    # If already looks absolute, leave it; else join with root
    if (grepl("^[A-Za-z]:[/\\\\]", x) || substr(x, 1, 1) == "/") {
      add_trailing_slash(x)
    } else {
      add_trailing_slash(file.path(root, x))
    }
  } else {
    x
  }
}

dir.list <- make_abs_paths(cfg$dir, root.path)

# Recursively create directories
create_dirs <- function(x) {
  if (is.list(x)) {
    invisible(lapply(x, create_dirs))
  } else if (is.character(x) && length(x) == 1L) {
    if (!dir.exists(x)) dir.create(x, recursive = TRUE)
  }
}

create_dirs(dir.list)

# Colors helper
get_colors <- function() {
  if (is.null(cfg$colors)) {
    stop("No 'colors' field found in config.yml")
  }
  cfg$colors
}

# Install/load comparatome
if (!requireNamespace("comparatome", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  devtools::install(paste0(tools.path, "comparatome"))
}
library(comparatome)

# GEO download settings from YAML
samples.opossum <- cfg$geo$samples_opossum
file.types      <- cfg$geo$file_types

samples.path <- file.path(dir.list$central, "samples")

for (sm in names(samples.opossum)) {
  sm.path <- file.path(samples.path, sm)
  if (!dir.exists(sm.path)) dir.create(sm.path, recursive = TRUE)
  
  geo.stem <- paste0(samples.opossum[[sm]], "_snRNA-seq_", sm)
  
  for (ft in file.types) {
    fname <- paste0(geo.stem, "_", ft)
    fpath <- file.path(sm.path, ft)
    url <- paste0(
      "https://www.ncbi.nlm.nih.gov/geo/download/?acc=",
      samples.opossum[[sm]],
      "&format=file&file=",
      URLencodeNCBI(fname)    # from comparatome
    )
    # system() call only if you actually want to auto-download on this machine
    # system(paste0("wget -O ", shQuote(fpath), " ", shQuote(url)))
  }
}
