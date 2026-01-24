# Compress large FASTA / GTF files for inst/extdata
# Run once during package development

input_dir <- "inst/data"
output_dir <- "inst/extdata"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(
    input_dir,
    pattern = "\\.(fa|fasta|gtf)$",
    full.names = TRUE
)

for (f in files) {
    out <- file.path(output_dir, paste0(basename(f), ".gz"))

    if (file.exists(out)) {
        message("Skipping existing: ", basename(out))
        next
    }

    message("Compressing: ", basename(f))
    con_in <- file(f, open = "rb")
    con_out <- gzfile(out, open = "wb")

    repeat {
        bytes <- readBin(con_in, what = raw(), n = 1e6)
        if (length(bytes) == 0) break
        writeBin(bytes, con_out)
    }

    close(con_in)
    close(con_out)
}

message("Compression complete.")
