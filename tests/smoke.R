library(dplyr)
dir.create("tests/tmp", showWarnings=FALSE)
input <- normalizePath("tests/tmp", winslash="/")
writexl::write_xlsx(list(VEPs=data.frame(
  variant_id_hg38=paste0("synthetic_",1:20),
  Type=rep(c("Pathogenic","Control"),each=10),
  SIFT_score=c(seq(.001,.045,length.out=10),seq(.1,.9,length.out=10)),
  REVEL_score=c(seq(.6,.95,length.out=10),seq(.01,.4,length.out=10))
)), file.path(input,"synthetic.xlsx"))
rscript <- file.path(R.home("bin"),"Rscript")
for (entry in c("ROC.R","Confusion Matrix.R","Heatmap.R")) {
  status <- system2(rscript,c(shQuote(entry),
    shQuote(paste0("--input=",file.path(input,"synthetic.xlsx"))),
    "--gene=synthetic", "--bootstrap=10"))
  stopifnot(status == 0)
  out <- "synthetic_benchmark_v5_results"
  workbook <- file.path(out,"01_tables/synthetic_benchmark_v5_summary.xlsx")
  stopifnot(file.exists(workbook))
  sheets <- readxl::excel_sheets(workbook)
  stopifnot(all(c("ROCSummary","BootstrapValidation","ThresholdMetrics","DbnsfpLatestVsBenchmark") %in% sheets))
  metrics <- readxl::read_excel(workbook,sheet="ThresholdMetrics")
  stopifnot(nrow(metrics)==2, all(metrics$TP==10),all(metrics$TN==10))
  roc <- readxl::read_excel(workbook,sheet="ROCSummary")
  stopifnot(all(roc$AUC==1), all(roc$AUPRC==1))
  stopifnot(length(list.files(out,pattern="\\.svg$",recursive=TRUE)) >= 5)
}
unlink("tests/tmp", recursive=TRUE)
unlink("synthetic_benchmark_v5_results",recursive=TRUE)
message("PASS: all three entries, workbook metrics, bootstrap, reference comparison and SVG outputs")
