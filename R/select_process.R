#' Select processing function to call
#'
#' Simple function to convert the method name (character) to the right function
#'
#' @param method character. The name of the method for DAF
#' @return a function performing the chosen analysis

select_process <- function(method){

    method <- tolower(method)

    f <- switch(method,
           deseq2 = process_DESeq2,
           metagenomeseq = process_metagenomeSeq,
           aldex2 = process_ALDEx2,
           mbzinb = process_mbzinb,
           edger = process_edgeR,
           raida = process_RAIDA,
           zibseq = process_ZIBseq,
           voom = process_voom,
           wilcox = process_wilcox)

    return(f)
}
