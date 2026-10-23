# Details ---------------------------------------------------------------
#      AUTHOR:	James Foster              DATE: 2026 10 23
#
# DESCRIPTION:  Export the simulated datasets from "Circular Dataset Structures.Rmd"
#               as CSV files, so that the same data can be modelled in Matlab
#               (../Matlab/CircularDatasetStructures.m).
#
#      INPUTS:  Circular Dataset Structures.Rmd, unwrap_functions.R
#
#     OUTPUTS:  ../Matlab/circular_datasets/*.csv
#
#	   CHANGES: -
#
#   REFERENCES: -
#
#    EXAMPLES:  Rscript "Export Circular Datasets.R" (run from the R folder)
#

# Extract the Rmd code ----------------------------------------------------
#only the chunks that simulate data are run in full,
#for the model fitting chunks only the data.frame assignments (dt_...) are run
rmd_file = 'Circular Dataset Structures.Rmd'
out_dir = file.path('..', 'Matlab', 'circular_datasets')
dir.create(out_dir, showWarnings = FALSE)

r_file = tempfile(fileext = '.R')
invisible(knitr::purl(input = rmd_file,
            output = r_file,
            documentation = 1,
            quiet = TRUE))
r_lines = readLines(r_file, encoding = 'UTF-8')
#chunk headers look like "## ----label-----..."
hdr = grep(pattern = '^## ----', x = r_lines)
chunk_labels = trimws(gsub(pattern = '^## -+|-+$', replacement = '', x = r_lines[hdr]))
chunk_code = mapply(FUN = function(from, to){r_lines[from:to]},
                    from = hdr + 1,
                    to = c(hdr[-1] - 1, length(r_lines)),
                    SIMPLIFY = FALSE)
names(chunk_code) = chunk_labels

#chunks to run in full (in document order)
full_chunks = c('Load functions and packages',
                'Set up plot colours',
                'Generate diverging data',
                'Simulate a change in direction',
                'Simulate a decrease in concentration',
                'Simulate high interindividual correlation',
                'Simulate data with low inter-individual heading correlation',
                'Simulate individuals with variable parameters',
                'Simulate individual heading changes')
stopifnot(all(full_chunks %in% chunk_labels))

# Run the simulations -----------------------------------------------------
pdf(file = NULL) # plots are not needed
for(lbl in chunk_labels)
{
  exprs = parse(text = chunk_code[[lbl]], encoding = 'UTF-8')
  if(lbl %in% full_chunks)
  {
    for(e in exprs){ eval(e, envir = globalenv()) }
  }else
  {
    #only assignments of datasets
    for(e in exprs)
    {
      if(is.call(e) &&
         (identical(e[[1]], as.name('=')) || identical(e[[1]], as.name('<-'))) &&
         grepl(pattern = '^dt_', x = deparse(e[[2]])) )
      {
        tryCatch(eval(e, envir = globalenv()),
                 error = function(err){}) # e.g. dt_aligned needs dt_divergence first
      }
    }
  }
}
invisible(dev.off())

# Save datasets -----------------------------------------------------------
dt_names = c('dt_divergence', 'dt_delta', 'dt_conc',
             'dt_highcorr', 'dt_lowcorr', 'dt_var', 'dt_hd')
for(nm in dt_names)
{
  write.csv(x = get(nm),
            file = file.path(out_dir, paste0(sub('^dt_', '', nm), '.csv')),
            row.names = FALSE)
}

#true parameters for individuals, for comparison with model estimates (degrees)
write.csv(x = data.frame(ID = seq_along(mu_highcorr), mu = as.numeric(mu_highcorr)),
          file = file.path(out_dir, 'highcorr_true.csv'), row.names = FALSE)
write.csv(x = data.frame(ID = seq_along(mu_lowcorr), mu = as.numeric(mu_lowcorr)),
          file = file.path(out_dir, 'lowcorr_true.csv'), row.names = FALSE)
write.csv(x = data.frame(ID = seq_along(mu_var), mu = as.numeric(mu_var),
                         kappa = kappa_id_var),
          file = file.path(out_dir, 'var_true.csv'), row.names = FALSE)
write.csv(x = data.frame(ID = seq_along(mu_hd), mu = as.numeric(mu_hd),
                         mu_delta = as.numeric(mu_delta),
                         kappa = kappa_id_hd),
          file = file.path(out_dir, 'hd_true.csv'), row.names = FALSE)

message('Datasets saved to ', normalizePath(out_dir))
