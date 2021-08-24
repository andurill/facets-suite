#' Get gene-level changes in copy number from FACETS output.
#'
#' This maps protein-coding genes onto copy-number segmentation from FACETS output derives the copy-number status from that, based on the integer copy number.
#' Additionally, the copy-number log-ratio of the SNPs falling insided the gene boundaries are used to perform a Z-test against the dipLogR baseline.
#' 
#' @param hisens_output Full FACETS hisensitivity output from \code{run_facets}.
#' @param purity_output Full FACETS purity output from \code{run_facets}.
#' @param gene_level    Full gene level output from \code{gene-level-changes}.
#' @param arm_level     Full arm level output from \code{arm-level-changes}.
#' @return \code{json}
#' 
#' @import data.table
#' @import jsonlite
#' 
#' @examples
#' \dontrun{
#' generate-json-ouput(hisens_ouput, NULL, gene_level_ouput, arm_level_output)
#' }

#' @export
generate_json = function(hisens_output,
                            purity_output,
                            gene_level,
                            arm_level) {
    saveRDS(list(hisens_output, purity_output, gene_level, arm_level), file="facets_test.Rdata")
    jsonlite::write_json(
        list("FACETS_PROCESSED_DATA_FOR_PLOT" = list(
            "VERTICAL_LINE_X_COORD_LIST" = "VERTICAL_LINE_X_COORD_LIST", 
            "Tick_Values" = "Tick_Values")), pretty=T,
        path = "json_output.txt")
    return()
    

    # grch37_coordinate 
    VERTICAL_LINE_X_COORD_LIST = c(0, 249250621, 492449994, 690472424, 881626700, 
                                    1062541960, 1233657027, 1392795690, 1539159712, 
                                    1680373143, 1815907890, 1950914406, 2084766301, 
                                    2199936179, 2307285719, 2409817111, 2500171864, 
                                    2581367074, 2659444322, 2718573305, 2781598825, 
                                    2829728720, 2881033286, 3036303846)

    # tick values
    Tick_Values = floor(head(VERTICAL_LINE_X_COORD_LIST, -1) + 
                        (tail(VERTICAL_LINE_X_COORD_LIST, -1) - 
                        head(VERTICAL_LINE_X_COORD_LIST, -1))/2)

    hisens = as.data.table(hisens_output[[1]]$segs)
    Line_Data = list()
    Segment_mean = list()
    for (i in seq(nrow(hisens))) {
        Line_Data[[i]] = list(
            "id" = unbox(paste0("Chromosome ", hisens[i,]$chrom, ": Copy #", i)),
            "data" = list(list("x" = unbox(hisens[i,]$start), "y" = unbox(hisens[i,]$cnlr.median)),
                    list("x" = unbox(hisens[i,]$end), "y" = unbox(hisens[i,]$cnlr.median))))
        
        seg_key = paste0("Chromosome", as.character(hisens[i,]$chrom), "Copy", i)
        Segment_mean[[i]] = unbox(hisens[i,]$cnlr.median)
        names(Segment_mean)[i] = paste0("Chromosome ", as.character(hisens[i,]$chrom), ": Copy #", i)
    }

    snps = as.data.table(hisens_output[[1]]$snps)
    Scatter_Data = list()
    for (i in seq(nrow(snps))) {
        Scatter_Data[[i]] = list(
            "id" = unbox(paste0("Chromosome ", hisens[i,]$chrom)),
            "data" = list(list("x" = unbox(hisens[i,]$maploc), "y" = unbox(hisens[i,]$cnlr.median)),
                        list("x" = unbox(hisens[i,]$maploc), "y" = unbox(hisens[i,]$cnlr.median))))
    }

    Line_Data_VALOR = list()
    Log_DR = c()
    Log_DR_counter = 1
    for (i in seq(nrow(hisens)-34)) {
        Line_Data_VALOR[[Log_DR_counter]] = list(
            "id" = unbox(paste0("Chromosome ", hisens[i,]$chrom, ": Copy #", i, "_Min")),
            "data" = list(list("x" = unbox(hisens[i,]$start), "y" = unbox(-sqrt(abs(hisens[i,]$mafR)))),
                            list("x" = unbox(hisens[i,]$end), "y" = unbox(-sqrt(abs(hisens[i,]$mafR))))))
        Line_Data_VALOR[[Log_DR_counter+1]] = list(
            "id" = unbox(paste0("Chromosome ", hisens[i,]$chrom, ": Copy #", i, "_Max")),
            "data" = list(list("x" = unbox(hisens[i,]$start), "y" = unbox(sqrt(abs(hisens[i,]$mafR)))),
                            list("x" = unbox(hisens[i,]$end), "y" = unbox(sqrt(abs(hisens[i,]$mafR))))))

        #seg_key_min = paste0("Chromosome ", as.character(hisens[i,]$chrom), "Copy #", i, "_Min")
        #seg_key_max = paste0("Chromosome ", as.character(hisens[i,]$chrom), "Copy #", i, "_Max")
        Log_DR[[Log_DR_counter]] = unbox(-sqrt(abs(hisens[i,]$mafR)))
        Log_DR[[Log_DR_counter+1]] = unbox(sqrt(abs(hisens[i,]$mafR)))
        names(Log_DR[Log_DR_counter]) = unbox(paste0("Chromosome ", as.character(hisens[i,]$chrom), "Copy #", i, "_Min"))
        names(Log_DR[Log_DR_counter+1]) = unbox(paste0("Chromosome ", as.character(hisens[i,]$chrom), "Copy #", i, "_Max"))
        Log_DR_counter = Log_DR_counter + 2
    }

    gene_dict = list()
for (i in seq(nrow(hisens))) {
  gene_list = strsplit(hisens[i,]$gene_list, ",")[[1]]
  gene_dict[[i]] = gene_list
  names(gene_dict)[i] = paste0("Chromosome ", as.character(hisens[i,]$chrom), ": Copy #", i)
  gene_line_data_cnlor = gene_line_data_valor = list()
  gene_scatter_data_cnlor = gene_scatter_data_valor = list()
  gene_line_data_CNCF = gene_line_data_EM = list()
  gene_cnlr_counter = 1
  gene_valor_counter = 1
  scatter_data_counter = 1
  for (g in gene_list) {
    scatter_data_for_gene = snps_gene_mapping[snps_gene_mapping$gene == g,]
    gene_start = unbox(VERTICAL_LINE_X_COORD_LIST[hisens[i,]$chrom] + genes[gene == g,]$gene_start)
    gene_end = unbox(VERTICAL_LINE_X_COORD_LIST[hisens[i,]$chrom] + genes[gene == g,]$gene_end)
    gene_cnlr_median = unbox(hisens[i,]$cnlr.median)
    gene_valor_positive_y = unbox(sqrt(abs(hisens[i,]$mafR)))
    gene_valor_negative_y = unbox(-sqrt(abs(hisens[i,]$mafR)))
    gene_line_data_cnlor[[gene_cnlr_counter]] = list(
      "id" = unbox(g),
      "data" = list(
        list(
          "x" = gene_start, 
          "y" = gene_cnlr_median),
        list(
          "x" = gene_end, 
          "y" = gene_cnlr_median)))
    gene_line_data_valor[[gene_valor_counter]] = list(
      "id" = unbox(paste0(g, " Copy #", as.character(gene_valor_counter), "_Min")),
      "data" = list(
        list(
          "x" = gene_start, 
          "y" = gene_valor_negative_y),
        list(
          "x" = gene_end, 
          "y" = gene_valor_negative_y)))
    gene_line_data_valor[[gene_valor_counter]] = list(
      "id" = unbox(paste0(g, " Copy #", as.character(gene_valor_counter), "_Max")),
      "data" = list(
        list(
          "x" = gene_start, 
          "y" = gene_valor_positive_y),
        list(
          "x" = gene_end, 
          "y" = gene_valor_positive_y)))
    gene_line_data_CNCF[[gene_cnlr_counter]] = list(
      "id" = unbox(paste0(g, " Copy #", as.character(gene_cnlr_counter), "_TCN")),
      "data" = list(
        list(
          "x" = gene_start, 
          "y" = gene_tcn,
        list(
          "x" = gene_end, 
          "y" = gene_tcn))))
    gene_line_data_CNCF[[gene_cnlr_counter]] = list(
      "id" = unbox(paste0(g, " Copy #", as.character(gene_cnlr_counter), "_LCN")),
      "data" = list(
        list(
          "x" = gene_start, 
          "y" = gene_lcn,
          list(
            "x" = gene_end, 
            "y" = gene_lcn))))
    gene_cnlr_counter = gene_cnlr_counter + 1
    gene_valor_counter= gene_valor_counter + 2
  }
}

    jsonlite::write_json(
        list("FACETS_PROCESSED_DATA_FOR_PLOT" = list(
            "VERTICAL_LINE_X_COORD_LIST" = VERTICAL_LINE_X_COORD_LIST, 
            "Tick_Values" = Tick_Values), 
            "CNLRPlot" = list("Line_Data" =  Line_Data_CNLOR,
                            #"Segment_mean" = Segment_mean,
                            "Scatter_Data" = Scatter_Data),
            "VALORPlot" = list("Line_Data" = Line_Data_VALOR,
                            "Log_DR" = Log_DR)), pretty=T,
        path = "json_output.txt")
}



