library(rentrez)
library(pbapply)

# very important
httr::set_config(httr::config(http_version = 0))

download_sequences <- function(ids, result_path = "results") {
  pblapply(ids,  function(ith_id) {
    # searching IDs of a specific conting
    search_res <- entrez_search(db = "nucleotide" , term = paste0(ith_id, "[Primary Accession]"))[["ids"]]
    
    # seaching the assembly to which the conting belongs
    assembly_link <- entrez_link(dbfrom = "nucleotide", db = "assembly" , id = search_res)[["links"]][["nuccore_assembly"]]
    
    # assembly accesion number
    assembly_summary <- entrez_summary(db = "assembly", id = assembly_link)[["assemblyaccession"]]
    
    # searching for all sequences belonging to the assembly
    all_seqs_wh <- entrez_search(db = "nucleotide" , term = paste0(assembly_summary, "[Assembly]"), 
                                 retmax = 0, use_history = TRUE)
    
    # downloading sequences
    all_seqs_fasta <- entrez_fetch(db = "nucleotide" , rettype = "fasta", web_history = all_seqs_wh[["web_history"]])
    
    # saving as fasta file
    cat(all_seqs_fasta, file = paste0(result_path, "/", assembly_summary, ".fa"))
    invisible(TRUE)
  })
}

download_sequences(c("NZ_UILC01000024.1", "NZ_UIJM01000020.1"))
