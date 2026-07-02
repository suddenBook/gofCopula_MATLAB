#!/usr/bin/env Rscript

# Export installed gofCopula datasets as raw IEEE 754 binary64 values.
# Usage: Rscript tools/exportRData.R OUTPUT_FOLDER

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: exportRData.R OUTPUT_FOLDER")
}

output_folder <- args[[1L]]
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

write_double <- function(value, destination) {
  connection <- file(destination, "wb")
  on.exit(close(connection))
  writeBin(as.double(value), connection, size = 8L, endian = "little")
}

for (dataset in c("Banks", "CryptoCurrencies")) {
  environment <- new.env(parent = emptyenv())
  data(list = dataset, package = "gofCopula", envir = environment)
  value <- get(dataset, envir = environment)

  for (year in names(value)) {
    write_double(value[[year]], file.path(
      output_folder, paste0(dataset, "_", year, ".bin")))
  }
  writeLines(names(value), file.path(
    output_folder, paste0(dataset, "_years.txt")))
  writeLines(colnames(value[[1L]]), file.path(
    output_folder, paste0(dataset, "_variables.txt")))
  writeLines(as.character(vapply(value, nrow, integer(1L))), file.path(
    output_folder, paste0(dataset, "_rows.txt")))
}

for (dataset in c("IndexReturns2D", "IndexReturns3D")) {
  environment <- new.env(parent = emptyenv())
  data(list = dataset, package = "gofCopula", envir = environment)
  value <- get(dataset, envir = environment)

  write_double(value, file.path(output_folder, paste0(dataset, ".bin")))
  writeLines(colnames(value), file.path(
    output_folder, paste0(dataset, "_variables.txt")))
  writeLines(as.character(dim(value)), file.path(
    output_folder, paste0(dataset, "_dim.txt")))
}
