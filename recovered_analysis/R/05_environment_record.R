software_versions <- data.frame(
  Software = c("R", "Seurat", "SeuratObject", "ggplot2"),
  Version = c(as.character(getRversion()), as.character(packageVersion("Seurat")),
              as.character(packageVersion("SeuratObject")), as.character(packageVersion("ggplot2"))))
save_table(software_versions, "software_versions.csv")
writeLines(capture.output(sessionInfo()), file.path(run_dir, "logs", "session_info.txt"))
writeLines(c("Mode:", run_mode, "Completed at:", as.character(Sys.time()),
             "Recovered source files:", normalizePath(recovery_dir, winslash = "/")),
           file.path(run_dir, "logs", "run_record.txt"))
message("Completed. NEW output folder: ", normalizePath(run_dir, winslash = "/"))
