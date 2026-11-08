###############################################################################
# run_all.R -- reproduces every result, table and figure (run from code/)
###############################################################################
Sys.setenv(LANG = "C.UTF-8", OPENBLAS_NUM_THREADS = "1")  # single-core timings
for (s in c("01_prepare_real_data.R", "02_simulation.R", "03_real_data.R", "04_scaling.R",
            "05_sensitivity.R", "06_equivalence_checks.R", "08_added_methods.R", "09_mingrm_diagnostics.R", "12_yq_exact.R", "13_yq_fe_toy.R", "14_yq_fe_all.R", "15_yq_fe_scatter.R", "16_information_split.R", "17_yq_fen_all.R", "18_yq_fen_scaling.R", "19_yq_fen_summary.R",
            # v12 (revision): evaluation by scenario and validation design, LODO, benchmark, aggregation, summaries
            "20_revision_eval.R", "22_distributed_aggregation.R", "21_timing_benchmark.R",
            "24_real_data_description.R", "23_revision_summary.R", "07_figures_tables.R",
            # v13: YQ-FE-Q (Supplementary Section S10)
            "25_yq_fe_q.R")) {
  message("==== ", s); t0 <- proc.time()[["elapsed"]]
  if (s %in% c("14_yq_fe_all.R", "17_yq_fen_all.R")) {   # two parts: simulation and real data
    for (part in c("sim", "real")) system2("Rscript", c(s, part))
  } else if (s == "25_yq_fe_q.R") {                       # six parts
    for (part in c("toy", "kern", "sim", "real", "lowrank", "fig")) system2("Rscript", c(s, part))
  } else if (s == "20_revision_eval.R") {                 # four parts (resumable)
    for (part in c("orig_sim", "orig_real", "new_sim", "lodo")) system2("Rscript", c(s, part))
  } else source(s, local = new.env())
  message(sprintf("==== %s finished in %.0f s", s, proc.time()[["elapsed"]] - t0))
}
writeLines(capture.output(sessionInfo()), file.path("..", "outputs", "logs", "sessionInfo_v12.txt"))
