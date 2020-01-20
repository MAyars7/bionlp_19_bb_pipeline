library(easyPubMed)

#config file
#Contains: my_dest_dir_prefix - output destination directory prefix
#          query_terms - summary of query terms. Ex: bacteria_habitat.
#          query_notes - Notes on query string generation.
#          site_hits - Number of hits returned by searching query on Pubmed site.
#          my_query_string - Pubmed-generated string to search for.
#          batch_size - Number of entries to store in each output file.
#          format - Format to write output in (XML, medline, abstract, etc.)

args <- commandArgs(trailingOnly = TRUE)
source(args[1])

#date for file labels
today <- Sys.Date()
dir_date <- format(today, format="%y%m%d")

my_dest_dir = sprintf("%s/%s_%s", my_dest_dir_prefix, query_terms, dir_date)
log_file = sprintf("%s/%s_%s.log", getwd(), query_terms, dir_date)

#Logging
sink(log_file, split=TRUE)

write(sprintf("Output destination directory: %s", my_dest_dir), stdout())
write(sprintf("Log file path: %s", log_file), stdout())
write(sprintf("Query notes: %s", query_notes), stdout())
write(sprintf("Query hits on Pubmed site: %s", site_hits), stdout())

dir.create(my_dest_dir)
out <- batch_pubmed_download(pubmed_query_string = my_query_string, batch_size=batch_size, format=output_format, dest_dir = my_dest_dir)

sink()