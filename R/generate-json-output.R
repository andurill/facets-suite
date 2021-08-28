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
    VERTICAL_LINE_X_COORD_LIST = c(
      0, 249250621, 492449994, 690472424, 881626700, 
      1062541960, 1233657027, 1392795690, 1539159712, 
      1680373143, 1815907890, 1950914406, 2084766301, 
      2199936179, 2307285719, 2409817111, 2500171864, 
      2581367074, 2659444322, 2718573305, 2781598825, 
      2829728720, 2881033286, 3036303846)

    # x-axis tick values for chromosomes
    Tick_Values = floor(head(VERTICAL_LINE_X_COORD_LIST, -1) +
                          (tail(VERTICAL_LINE_X_COORD_LIST, -1) - 
                            head(VERTICAL_LINE_X_COORD_LIST, -1))/2)


    # seg data for line plots
    segs = as.data.table(hisens_output[[1]]$segs)

    # gene data for line plots
    genes = as.data.table(gene_level[[1]])

    # snp level data for scatter plots
    snps = as.data.table(hisens_output[[1]]$snps)
    # assign dummy end location, which will be
    #  same as start location. Both start and end 
    #  are required for foverlaps function
    snps$maploc_end = snps$maploc
    # assing snp to a gene base on genomics coordinates
    setkey(genes, chrom, gene_start, gene_end)
    setkey(snps, chrom, maploc, maploc_end)
    snps_gene_mapping = data.table::foverlaps(snps, genes, type="within")

    # seg data for line plots mapped to all the genes
    #  within a given segment
    hisens = as.data.table(hisens_output[[1]]$segs)
    hisens[aggregate(genes[,c(2)], 
                    by = genes[,c(3,10,11)], function(genes) paste(genes, collapse = ",")), 
          on=c(chrom = "chrom", 
                start = "seg_start", 
                end = "seg_end"), gene_list := i.gene]

    #hisens_by_chr = hisens[, list(start = min(start), end = max(end)), by = chrom]


    # Log odds plot line data
    Line_Data_CNLOR = list()
    Segment_mean = list()
    for (i in seq(nrow(hisens))) {
      Line_Data_CNLOR[[i]] = list(
        "id" = unbox(paste0("Chromosome ", hisens[i,]$chrom, ": Copy #", i)),
        "data" = list(
          list(
            "x" = unbox(VERTICAL_LINE_X_COORD_LIST[hisens[i,]$chrom] + hisens[i,]$start), 
            "y" = unbox(hisens[i,]$cnlr.median)),
          list(
            "x" = unbox(VERTICAL_LINE_X_COORD_LIST[hisens[i,]$chrom] + hisens[i,]$end), 
            "y" = unbox(hisens[i,]$cnlr.median))))
      
      #seg_key = paste0("Chromosome ", as.character(hisens[i,]$chrom), " Copy #", i)
      Segment_mean[[i]] = unbox(hisens[i,]$cnlr.median)
      names(Segment_mean)[i] = paste0("Chromosome ", as.character(hisens[i,]$chrom), ": Copy #", i)
    }

    # Log odds and variant log odds plot scatter data
    snps = as.data.table(hisens_output[[1]]$snps)
    CNLR_Scatter_Data = VALOR_Scatter_Data = list()
    for (chrom_var in unique(snps$chrom)) {
      cnlr_data_points = valor_data_points = list()
      snps_per_chrom = snps[snps$chrom == chrom_var,]
      for (i in seq(1, nrow(snps_per_chrom))) {
        cnlr_data_points[[i]] = list(
          "x" = unbox(VERTICAL_LINE_X_COORD_LIST[chrom_var] + snps_per_chrom[i,]$maploc), 
          "y" = unbox(snps_per_chrom[i,]$cnlr))
        if (!is.na(snps_per_chrom[i,]$valor)) {
          valor_data_points[[i]] = list(
            "x" = unbox(VERTICAL_LINE_X_COORD_LIST[chrom_var] + snps_per_chrom[i,]$maploc), 
            "y" = unbox(snps_per_chrom[i,]$valor))
        }
        CNLR_Scatter_Data[[chrom_var]] = list(
          "id" = unbox(paste0("Chromosome ", chrom_var)),
          "data" = cnlr_data_points
        )
        VALOR_Scatter_Data[[chrom_var]] = list(
          "id" = unbox(paste0("Chromosome ", chrom_var)),
          "data" = valor_data_points
        )
      }
    }

    # variant log odds ratio plot line data
    Line_Data_VALOR = Log_DR = list()
    Log_DR_counter = 1
    for (i in seq(nrow(hisens))) {
      Line_Data_VALOR[[Log_DR_counter]] = list(
          "id" = unbox(paste0("Chromosome ", hisens[i,]$chrom, ": Copy #", i, "_Min")),
          "data" = list(
            list(
              "x" = unbox(VERTICAL_LINE_X_COORD_LIST[hisens[i,]$chrom] + hisens[i,]$start), 
              "y" = unbox(-sqrt(abs(hisens[i,]$mafR)))),
            list(
              "x" = unbox(VERTICAL_LINE_X_COORD_LIST[hisens[i,]$chrom] + hisens[i,]$end), 
              "y" = unbox(-sqrt(abs(hisens[i,]$mafR))))))
      Line_Data_VALOR[[Log_DR_counter+1]] = list(
          "id" = unbox(paste0("Chromosome ", hisens[i,]$chrom, ": Copy #", i, "_Max")),
          "data" = list(
            list(
              "x" = unbox(VERTICAL_LINE_X_COORD_LIST[hisens[i,]$chrom] + hisens[i,]$start), 
              "y" = unbox(sqrt(abs(hisens[i,]$mafR)))),
            list(
              "x" = unbox(VERTICAL_LINE_X_COORD_LIST[hisens[i,]$chrom] + hisens[i,]$end), 
              "y" = unbox(sqrt(abs(hisens[i,]$mafR))))))

      Log_DR[[Log_DR_counter]] = unbox(-sqrt(abs(hisens[i,]$mafR)))
      Log_DR[[Log_DR_counter+1]] = unbox(sqrt(abs(hisens[i,]$mafR)))
      names(Log_DR)[Log_DR_counter] = paste0("Chromosome ", as.character(hisens[i,]$chrom), ": Copy #", i, "_Min")
      names(Log_DR)[Log_DR_counter+1] = paste0("Chromosome ", as.character(hisens[i,]$chrom), ": Copy #", i, "_Max")
      Log_DR_counter = Log_DR_counter + 2
    }

    # integer copy number line data
    CNCF_ICNP = EM_ICNP = list()
    ICNP_counter = 1
    for (i in seq(nrow(segs))) {
      CNCF_ICNP[[ICNP_counter]] = list(
        "id" = unbox(paste0("Chromosome ", segs[i,]$chrom, ": Copy #", i, "_TCN")),
        "data" = list(
          list(
            "x" = unbox(VERTICAL_LINE_X_COORD_LIST[segs[i,]$chrom] + segs[i,]$start), 
            "y" = unbox(segs[i,]$tcn)),
          list(
            "x" = unbox(VERTICAL_LINE_X_COORD_LIST[segs[i,]$chrom] + segs[i,]$end), 
            "y" = unbox(segs[i,]$tcn))))
      CNCF_ICNP[[ICNP_counter+1]] = list(
        "id" = unbox(paste0("Chromosome ", segs[i,]$chrom, ": Copy #", i, "_LCN")),
        "data" = list(
          list(
            "x" = unbox(VERTICAL_LINE_X_COORD_LIST[segs[i,]$chrom] + segs[i,]$start), 
            "y" = unbox(segs[i,]$lcn)),
          list(
            "x" = unbox(VERTICAL_LINE_X_COORD_LIST[segs[i,]$chrom] + segs[i,]$end), 
            "y" = unbox(segs[i,]$lcn))))
      EM_ICNP[[ICNP_counter]] = list(
        "id" = unbox(paste0("Chromosome ", segs[i,]$chrom, ": Copy #", i, "_TCN_EM")),
        "data" = list(
          list(
            "x" = unbox(VERTICAL_LINE_X_COORD_LIST[segs[i,]$chrom] + segs[i,]$start), 
            "y" = unbox(segs[i,]$tcn.em)),
          list(
            "x" = unbox(VERTICAL_LINE_X_COORD_LIST[segs[i,]$chrom] + segs[i,]$end), 
            "y" = unbox(segs[i,]$tcn.em))))
      EM_ICNP[[ICNP_counter+1]] = list(
        "id" = unbox(paste0("Chromosome ", segs[i,]$chrom, ": Copy #", i, "_LCN_EM")),
        "data" = list(
          list(
            "x" = unbox(VERTICAL_LINE_X_COORD_LIST[segs[i,]$chrom] + segs[i,]$start), 
            "y" = unbox(segs[i,]$lcn.em)),
          list(
            "x" = unbox(VERTICAL_LINE_X_COORD_LIST[segs[i,]$chrom] + segs[i,]$end), 
            "y" = unbox(segs[i,]$lcn.em))))
      Log_DR_counter = Log_DR_counter + 2
    }
      

    gene_dict = list()
    gene_dict_counter = 1
    for (i in seq(nrow(hisens))) {
      gene_list = strsplit(hisens[i,]$gene_list, ",")[[1]]
      # a segment has no gene, skip and go to the next
      if(all(is.na(gene_list))) {
        next
      }
      gene_data = list()
      gene_counter = 1
      for (g in gene_list) {
        scatter_data_for_gene = snps_gene_mapping[snps_gene_mapping$gene == g,]
        # snp level data organization for each gene
        # scatter_data_for_gene should never be empty for a given segment
        scatter_data_counter = 1
        for (s in seq(1, nrows(scatter_data_for_gene))) {
          gene_scatter_data_cnlor[[scatter_data_counter]] = list(
            "x" = unbox(VERTICAL_LINE_X_COORD_LIST[hisens[i,]$chrom] + scatter_data_for_gene[s,]$maploc), 
            "y" = unbox(scatter_data_for_gene[s,]$cnlr))
          if (!is.na(scatter_data_for_gene[s,]$valor)) {
            gene_scatter_data_valor[[scatter_data_counter]] = list(
              "x" = unbox(VERTICAL_LINE_X_COORD_LIST[hisens[i,]$chrom] + scatter_data_for_gene[s,]$maploc), 
              "y" = unbox(scatter_data_for_gene[s,]$valor))
          }
          scatter_data_counter = scatter_data_counter + 1
        }
        # gene level data organization
        gene_start = unbox(VERTICAL_LINE_X_COORD_LIST[hisens[i,]$chrom] + genes[gene == g,]$gene_start)
        gene_end = unbox(VERTICAL_LINE_X_COORD_LIST[hisens[i,]$chrom] + genes[gene == g,]$gene_end)
        gene_tcn = unbox(hisens[i,]$tcn)
        gene_lcn = unbox(hisens[i,]$lcn)
        gene_tcn.em = unbox(hisens[i,]$tcn.em)
        gene_lcn.em = unbox(hisens[i,]$lcn.em)
        gene_cnlr_median = unbox(hisens[i,]$cnlr.median)
        gene_valor_positive_y = unbox(sqrt(abs(hisens[i,]$mafR)))
        gene_valor_negative_y = unbox(-sqrt(abs(hisens[i,]$mafR)))
        gene_line_data_cnlor = list(
          "id" = unbox(g),
          "data" = list(
            list(
              "x" = gene_start, 
              "y" = gene_cnlr_median),
            list(
              "x" = gene_end, 
              "y" = gene_cnlr_median)))
        gene_line_data_valor_min = list(
          "id" = unbox(paste0(g, " Copy #", as.character(gene_valor_counter), "_Min")),
          "data" = list(
            list(
              "x" = gene_start, 
              "y" = gene_valor_negative_y),
            list(
              "x" = gene_end, 
              "y" = gene_valor_negative_y)))
        gene_line_data_valor_max = list(
          "id" = unbox(paste0(g, " Copy #", as.character(gene_valor_counter), "_Max")),
          "data" = list(
            list(
              "x" = gene_start, 
              "y" = gene_valor_positive_y),
            list(
              "x" = gene_end, 
              "y" = gene_valor_positive_y)))
        gene_line_data_CNCF_TCN = list(
          "id" = unbox(paste0(g, " Copy #", as.character(gene_cnlr_counter), "_TCN")),
          "data" = list(
            list(
              "x" = gene_start, 
              "y" = gene_tcn,
              list(
                "x" = gene_end, 
                "y" = gene_tcn))))
        gene_line_data_CNCF_LCN = list(
          "id" = unbox(paste0(g, " Copy #", as.character(gene_cnlr_counter), "_LCN")),
          "data" = list(
            list(
              "x" = gene_start, 
              "y" = gene_lcn,
              list(
                "x" = gene_end, 
                "y" = gene_lcn))))
        gene_line_data_EM_TCN = list(
          "id" = unbox(paste0(g, " Copy #", as.character(gene_cnlr_counter), "_TCN_EM")),
          "data" = list(
            list(
              "x" = gene_start, 
              "y" = gene_tcn.em,
              list(
                "x" = gene_end, 
                "y" = gene_tcn.em))))
        gene_line_data_EM_LCN = list(
          "id" = unbox(paste0(g, " Copy #", as.character(gene_cnlr_counter), "_LCN_EM")),
          "data" = list(
            list(
              "x" = gene_start, 
              "y" = gene_lcn.em,
              list(
                "x" = gene_end, 
                "y" = gene_lcn.em))))
        
        # organize line and scatter data for the gene, g
        gene_data_current_gene = list(
          g = list(
            "CNLRPlot" = list(
              "Line_Data" = gene_line_data_cnlor,
              "Scatter_Data" = list(
                "id" = g,
                "data" = gene_scatter_data_cnlor)
            ),
            "VALORPlot" = list(
              "Line_Data" = list(
                gene_line_data_valor_min, 
                gene_line_data_valor_max),
              "Scatter_Data" = list(
                "id" = g,
                "data" = gene_scatter_data_valor)
            ),
            "ICNPlot" = list(
              "Line_Data" = list(
                "CNCF" = list(
                  gene_line_data_CNCF_TCN,
                  gene_line_data_CNCF_LCN),
                "EM" = list(
                  gene_line_data_EM_TCN,
                  gene_line_data_EM_LCN)
              ),
              "Scatter_Data" = NULL
            ),
            "Gene_Info" = unbox(scatter_data_for_gene$cn_state[1])
          )
        )
        gene_data[[gene_counter]] = gene_data_current_gene
        names(gene_data)[gene_counter] = unbox(g)
        gene_counter = gene_counter + 1
      }
      # organize data for all genes in a given segment
      segment_id = paste0("Chromosome ", as.character(hisens[i,]$chrom), ": Copy #", gene_dict_counter)
      gene_dict[[gene_dict_counter]] = list(
        "Gene_List" = gene_list,
        "Gene_Data" = gene_data)
      names(gene_dict)[gene_dict_counter] = unbox(segment_id)
    }

    jsonlite::write_json(
      list(
        "FACETS_PROCESSED_DATA_FOR_PLOT" = list(
          "VERTICAL_LINE_X_COORD_LIST" = VERTICAL_LINE_X_COORD_LIST, 
          "Tick_Values" = Tick_Values, 
          "CNLRPlot" = list("Line_Data" =  Line_Data_CNLOR,
                            "Segment_mean" = Segment_mean,
                            "Scatter_Data" = CNLR_Scatter_Data),
          "VALORPlot" = list("Line_Data" = Line_Data_VALOR,
                            "Log_DR" = Log_DR,
                            "Scatter_Data" = VALOR_Scatter_Data),
          "ICNPlot" = list("CNCF" = CNCF_ICNP,
                          "EM" = EM_ICNP),
          "Gene_Dictionary" = gene_dict,
          "All_Genes_List" = genes$gene),
        "Parameters" = list(
          "cval" = cval,
          "purity_cval" = purity_cval,
          "min_het" = min_het,
          "purity_min_het" = purity_min_het,
          "normal_depth" = normal_depth
        )
      ), pretty=T, path = "~/json.txt")
}



