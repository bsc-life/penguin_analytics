#libraries install and load
source("packages_installer.R")
# helper functions
source("functions.r")


ppi_filter = "filtered"
## ppi_filter: filtered , filtered Protein Proteins interactions to be LNCaP-specific and expressed in LNCaP cell line )
## ppi_filter: unfiltered , unfiltered Protein Proteins interactions include all
cell_line = "LNCaP" 
data_folder = file.path("src/EPIN_reconstruction_data", cell_line)
general_data_folder = file.path("src/EPIN_reconstruction_data")
output_folder = paste0("outputs_", cell_line, "_EPINS")
dir.create(output_folder,  showWarnings = FALSE)
dir.create(file.path(output_folder, "tables"),  showWarnings = FALSE)
dir.create(file.path(output_folder, "tables", ppi_filter),  showWarnings = FALSE)

#create output folders
for(subDir in c("EP_graph_edges", "TF_tables", "Intermediate_nodes_tables", "Promoter_genes_tables"))
{
  dir.create(file.path(output_folder, "tables", ppi_filter ,subDir), showWarnings = FALSE)
}



## variables
FPKM_threshold = 0.003 # both replicates to be above
sliding_window = 100  # sliding window for filtering fimo
FIMO_pvalue = 1e-4
max_num_edges = 2 # mux number of edges linking E to P in the PPI reconstruction



#####
gwas_genes <- as.data.frame(read.table(file.path(data_folder, "PrCa_GeneList_Used.csv"), sep = ",", header = T))

gwas_promoters <-  as.data.frame(read.table(file.path(data_folder, "GWAS_promoters.txt"), sep = "\t", header = T, stringsAsFactors = F))

SNPs <- as.data.frame(read.table( file.path(data_folder, "paintor_1causals.txt"), header = TRUE, stringsAsFactors = F)) %>%
  select(-c(BP, A0, A1, Z, BETA, region, index_rsid, n, SE, N_CONTROLS, N_CASES, ID, Posterior_Prob, ch_pos,N)) %>%
  rename(SNP_start = start, SNP_stop = stop)

##### READ the list of DNA binding proteins
motifs2gene <- as.data.frame(read.table(file.path(general_data_folder,  "JASPAR_motifs_2020_human.txt"), header = TRUE, sep = "\t", stringsAsFactors = F))

##### READ read counts 
expression <- as.data.frame(read.table(file.path(data_folder, "Cuff_Gene_Counts.txt"), header = T)) %>%
  filter(LNCaP_1 > FPKM_threshold & LNCaP_2 > FPKM_threshold)


#### READ PPIs. Downloaded PPIS from IID database
ppi_list <- as.data.frame(read.delim(file.path(general_data_folder, "human_annotated_PPIs_IID.txt"),  header = TRUE, stringsAsFactors = F)) %>% 
  select(uniprot1, uniprot2, symbol1, symbol2, prostate, prostate.carcinoma, prostate.cancer, nucleus, evidence.type, drug.targets, methods) %>%
  filter(grepl(pattern = "exp", evidence.type)) %>%
  filter(nucleus == 1 ) %>%
  mutate(num_methods = str_count(methods, ";") + 1) %>%
  filter(num_methods >= 2) %>%
  rename(to = symbol2, from = symbol1 )

if(ppi_filter == "filtered")
{
  ppi_list <- ppi_list %>%
    filter(prostate.carcinoma == 1 | prostate.cancer == 1 | prostate == 1 ) %>%
    filter(to %in% expression$Gene_ID & from %in% expression$Gene_ID ) %>%
    select(-methods) %>%
    select(from , to)
}else{
  ppi_list <- ppi_list %>%
    select(-methods) %>%
    select(from , to)
}


ppi_net <- simplify(graph_from_data_frame(ppi_list %>% select(from, to), directed=F))

#REEAD PROMOTER regions
promoters <- as.data.frame(read.table(file.path(data_folder, "genepos.txt"), sep = "\t", header = F)) %>%
  rename(seqnames = V1, start = V2, end = V3, gene = V4) %>%
  mutate(seqnames = paste("chr", seqnames, sep = ""))


## READ the refined enhancer anchors
loop_file = file.path(data_folder,  "E-P_loops_refined-regions.tsv")
loops <- as_tibble(read.table(loop_file, comment.char="", header = T,  sep = "\t", stringsAsFactors = F)) %>%
  mutate(enhancer_anchor_id = paste(chromosome.Enhancer.bin, start.Enhancer.bin, end.Enhancer.bin, sep = "_"),
         promoter_anchor_id = paste(chromosome.Promoter, start.Promoter, end.Promoter, sep = "_"),
         loopid = paste(enhancer_anchor_id, promoter_anchor_id, sep =  "_"),
         promoter_gene = gene.name) %>%
  distinct(gene.name, promoter_gene, loopid, enhancer_anchor_id, promoter_anchor_id, start.Enhancer.bin, chromosome.Enhancer.bin, end.Enhancer.bin, chromosome.Promoter , start.Promoter, end.Promoter) 



##### READ CTCF binding sites in LCNaP cell line
chipseq_CTCF <- as.data.frame(read.table(file.path(data_folder, "CTCF_lncap_ENCFF155SPQ.bed")))



make_graph_per_gene <- function(loop, max_num_edges, save_me = 1)
{
  ptm=proc.time()
  my_gene <- unique(loop$promoter_gene)

  ppi_sp <- c()
  p_dbp <- data.frame()
  e_dbp <- data.frame()
  
  
  my_enhancers <- loop %>% filter(gene.name == my_gene)
  if(dim(my_enhancers)[1] == 0){return(data.frame())}
  
  ##### PROMOTERS
  my_promoter <- promoters %>% filter(gene == my_gene)
  if(dim(my_promoter)[1] == 0){return(data.frame())}
  
  
  chromo_promo = str_replace(my_promoter$seqnames, "chr", "")
  p_fimo_file <- file.path(data_folder, "fimo_promoters/")
  p_fimo <- read_fimo(p_fimo_file, paste(chromo_promo,  my_promoter$start, my_promoter$end, sep = "_"), 0)
  
  ## add the chipseq_info of CTCF
  add_ctcf_chipseq_info <- function(fimo , chr , start, end)
  { 
    chipseq <- chipseq_CTCF %>% filter(V1 == chr &  V2 >= start &  V3 <= end) %>% 
      select(V1, V2, V3) %>% mutate(Chipseq_CTCF = 1) %>%
      rename(ctcf_start = V2, ctcf_chr = V1, ctcf_end = V3) %>% 
      as.data.table()
    
    if(dim(chipseq)[1] == 0){return(fimo %>% mutate(Chipseq_CTCF = NA))}
    setkey(chipseq, ctcf_chr, ctcf_start, ctcf_end)
    p_fimo <- as.data.frame(foverlaps(as.data.table(p_fimo), chipseq, by.x=c("chr_DNA", "start", "stop"), type="any")) %>%
      select(-ctcf_start,-ctcf_end )
    return(p_fimo)
  }
  
  p_fimo <- add_ctcf_chipseq_info(p_fimo, unique(my_promoter$seqnames), min(my_promoter$start),  min(my_promoter$end))

  p_dbp_current <- data.frame(TF = unique(p_fimo$motif_alt_id))
  p_dbp_current <- p_dbp_current %>% filter(TF %in% V(ppi_net)$name)
  if(dim(p_dbp_current)[1] == 0){return(data.frame())}
  p_dbp <- rbind(p_dbp, p_dbp_current)
  
  CTCF_ChipSeq_P = ifelse(dim(p_fimo %>% filter(!is.na(Chipseq_CTCF)))[1] > 0, 1 , 0)
  CTCF_P <- ifelse("CTCF" %in% p_fimo$motif_alt_id, 1 , 0)
  
  
  if(dim(p_fimo)[1] == 0){return(data.frame())}
  
  snp_in_P = data.frame(SNP_in_P = dim(p_fimo %>% filter( SNP_in_Motif == 1 ) %>% distinct(motif_alt_id, rsid) %>% mutate(anchor = "P"))[1])
  
  ##### ENHANCER(S)
  e_with_ctcf <- c()
  e_with_ctcf_chipseq <- c()

  
  for(enhancer in unique(my_enhancers$enhancer_anchor_id))
  {
    print(paste("enhancer", enhancer))
    
    e_fimo_file <- file.path(data_folder, "fimo_enhancers/")
    e_fimo <- read_fimo(e_fimo_file, enhancer, 0)
    
    if(dim(e_fimo)[1] == 0){next}
    
    e_fimo <- add_ctcf_chipseq_info(e_fimo, str_split(enhancer, "_")[[1]][1], str_split(enhancer, "_")[[1]][2], str_split(enhancer, "_")[[1]][3])
    

    e_dbp_current <- data.frame(TF = unique(e_fimo$motif_alt_id), position = enhancer)
    e_dbp_current <- e_dbp_current %>% filter(TF %in% V(ppi_net)$name)
    if(dim(e_dbp_current)[1] == 0){next}
    e_dbp <- rbind(e_dbp, e_dbp_current)
    
    CTCF_ChipSeq_E = ifelse(dim(e_fimo %>% filter(!is.na(Chipseq_CTCF)))[1] > 0, 1 , 0)
    CTCF_E <- ifelse("CTCF" %in% e_fimo$motif_alt_id, 1 , 0)
 
       ##
    if("CTCF" %in% e_fimo$motif_alt_id){
      e_with_ctcf <- c(e_with_ctcf, enhancer)}
    if(CTCF_ChipSeq_E > 0){
      e_with_ctcf_chipseq <- c(e_with_ctcf_chipseq, enhancer)}
    
    ppi_sp <- c()
    for(i in p_dbp_current$TF){ppi_sp <- c(ppi_sp, all_shortest_paths(ppi_net, from=i, to=as.character(e_dbp_current$TF))$res)}
    ppi_sp <- ppi_sp[lapply(ppi_sp, length) <= max_num_edges ]
    
    if(length(ppi_sp)  > 0 )
    {
      if(exists("net")) { net <- graph.union(net, find_network(ppi_sp, e_dbp_current, p_dbp)) }
      if(!exists("net")) { net <-  find_network(ppi_sp, e_dbp_current, p_dbp) }
    }
  }

  if(dim(e_dbp)[1] == 0){return(data.frame())}
  
  #### Delete an intermediate if it is just connecting two DBPs of the same anchor
  ValidIntermediate <- function(i){
    my_neighbors <- paste(names(i),collapse ="_" )
    if( str_count(my_neighbors, "ENHA") < length(i)  && str_count(my_neighbors, "PROM") < length(i))
    {
      return(TRUE) 
    }else{
      return(FALSE)
    }
  }
  intermediates <- V(net)$name[!grepl("ENHA_", V(net)$name) & !grepl("PROM_", V(net)$name) ]
  
  if(length(intermediates) > 0){
  
    neig <- adjacent_vertices(net, intermediates)
    valid_intermediates <- intermediates[which(sapply(neig, ValidIntermediate))]
    remove_intermediates <- intermediates[intermediates %nin% valid_intermediates ]
  
    net <- delete.vertices(net, remove_intermediates)
  
  }
  net <- delete.vertices(net, degree(net)==0)
  

  
  as_long_data_frame(net) %>% 
    select(3,4) %>%
    write.table(file.path(output_folder,  "tables", ppi_filter, "EP_graph_edges", paste(my_gene, "txt", sep = ".")), quote = F, row.names = F, sep = "\t", col.names = F)
  
  
  print(proc.time() - ptm)
  
}



for(gene in unique(loops$promoter_gene))
{
  print(gene)
  loop <- loops %>% filter(promoter_gene == gene)
  make_graph_per_gene(loop, max_num_edges)
}




