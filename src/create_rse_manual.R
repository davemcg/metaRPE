# bug fix for this function as there is a slight gene name 
# misalignment between the provided G029 GTF and the unify output
create_rse_manual2 <- function (project, project_home = project_homes(organism = organism, 
                                                                      recount3_url = recount3_url), type = c("gene", "exon", "jxn"), 
                                organism = c("human", "mouse"), annotation = annotation_options(organism), 
                                bfc = recount3_cache(), jxn_format = c("ALL", "UNIQUE"), 
                                recount3_url = getOption("recount3_url", "http://duffel.rail.bio/recount3"), 
                                verbose = getOption("recount3_verbose", TRUE)) 
{
  type <- match.arg(type)
  organism <- match.arg(organism)
  project_home <- match.arg(project_home)
  annotation <- match.arg(annotation)
  jxn_format <- match.arg(jxn_format)
  if (verbose) {
    message(Sys.time(), " downloading and reading the metadata.")
  }
  metadata <- read_metadata(file_retrieve(url = locate_url(project = project, 
                                                           project_home = project_home, type = "metadata", organism = organism, 
                                                           annotation = annotation, recount3_url = recount3_url), 
                                          bfc = bfc, verbose = verbose))
  project_home_original <- project_home
  project_home <- metadata$recount_project.file_source[1]
  metadata$BigWigURL <- locate_url(project = project, project_home = project_home, 
                                   type = "bw", organism = organism, annotation = annotation, 
                                   recount3_url = recount3_url, sample = metadata$external_id)
  if (type == "jxn") {
    jxn_files <- locate_url(project = project, project_home = project_home, 
                            type = "jxn", organism = organism, annotation = annotation, 
                            jxn_format = jxn_format, recount3_url = recount3_url)
  }
  if (verbose) {
    message(Sys.time(), " downloading and reading the feature information.")
  }
  if (type %in% c("gene", "exon")) {
    feature_info <- rtracklayer::import.gff(file_retrieve(url = locate_url_ann(type = type, 
                                                                               organism = organism, annotation = annotation, recount3_url = recount3_url), 
                                                          bfc = bfc, verbose = verbose))
  }
  else if (type == "jxn") {
    feature_info <- utils::read.delim(file_retrieve(url = jxn_files[grep("\\.RR\\.gz$", 
                                                                         jxn_files)], bfc = bfc, verbose = verbose))
    feature_info$strand[feature_info$strand == "?"] <- "*"
    feature_info <- GenomicRanges::GRanges(feature_info)
  }
  if (verbose) {
    message(Sys.time(), " downloading and reading the counts: ", 
            nrow(metadata), ifelse(nrow(metadata) > 1, " samples", 
                                   " sample"), " across ", length(feature_info), 
            " features.")
  }
  if (type %in% c("gene", "exon")) {
    counts <- read_counts(file_retrieve(url = locate_url(project = project, 
                                                         project_home = project_home, type = type, organism = organism, 
                                                         annotation = annotation, recount3_url = recount3_url), 
                                        bfc = bfc, verbose = verbose), samples = metadata$external_id)
  }
  else if (type == "jxn") {
    counts <- Matrix::readMM(file_retrieve(url = jxn_files[grep("\\.MM\\.gz$", 
                                                                jxn_files)], bfc = bfc, verbose = verbose))
    if (verbose) {
      message(Sys.time(), " matching exon-exon junction counts with the metadata.")
    }
    jxn_rail <- read.delim(file_retrieve(url = jxn_files[grep("\\.ID\\.gz$", 
                                                              jxn_files)], bfc = bfc, verbose = verbose))
    m <- match(metadata$rail_id, jxn_rail$rail_id)
    stopifnot(`Metadata rail_id and exon-exon junctions rail_id are not matching.` = !all(is.na(m)))
    counts <- counts[, m, drop = FALSE]
    colnames(counts) <- metadata$external_id
  }
  if (verbose) {
    message(Sys.time(), " construcing the RangedSummarizedExperiment (rse) object.")
  }
  stopifnot(`Metadata external_id and counts colnames are not matching.` = identical(metadata$external_id, 
                                                                                     colnames(counts)))
  if (type == "gene") {
    print("INVOKING DAVID VERSION")
    # reorder feature_info by rownames(counts) to align if possible
    names(feature_info) <- feature_info$gene_id
    feature_info <- feature_info[rownames(counts),]
    stopifnot(`Gene names and count rownames are not matching.` = identical(feature_info$gene_id, 
                                                                            rownames(counts)))
  }
  else if (type == "exon") {
    stopifnot(`Exon names and count rownames are not matching.` = identical(feature_info$recount_exon_id, 
                                                                            rownames(counts)))
  }
  else if (type == "jxn") {
    rownames(counts) <- as.character(feature_info)
  }
  names(feature_info) <- rownames(counts)
  rownames(metadata) <- colnames(counts)
  recount3_pkg <- sessioninfo::package_info(pkgs = "recount3", 
                                            include_base = FALSE, dependencies = FALSE)
  rse <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = counts), 
                                                    colData = S4Vectors::DataFrame(metadata, check.names = FALSE), 
                                                    rowRanges = feature_info, metadata = list(time_created = Sys.time(), 
                                                                                              recount3_version = as.data.frame(recount3_pkg), project = project, 
                                                                                              project_home = project_home_original, type = type, 
                                                                                              organism = organism, annotation = annotation, jxn_format = jxn_format, 
                                                                                              recount3_url = recount3_url))
  if (type %in% c("gene", "exon")) {
    assayNames(rse) <- "raw_counts"
    metadata(rse)$jxn_format <- NULL
  }
  return(rse)
}

# this is a direct copy expect I call my modified create_rse_manual2
create_rse2 <- function (project_info, type = c("gene", "exon", "jxn"), annotation = annotation_options(project_info$organism), 
                         bfc = recount3_cache(), jxn_format = c("ALL", "UNIQUE"), 
                         recount3_url = getOption("recount3_url", "http://duffel.rail.bio/recount3"), 
                         verbose = getOption("recount3_verbose", TRUE)) 
{
  stopifnot(`'project_info' should be a data.frame` = is.data.frame(project_info), 
            `'project_info' should only have one row` = nrow(project_info) == 
              1, `'project_info' should contain columns 'project', 'project_home' and 'organism'.` = all(c("project", 
                                                                                                           "project_home", "organism") %in% colnames(project_info)))
  type <- match.arg(type)
  annotation <- match.arg(annotation)
  jxn_format <- match.arg(jxn_format)
  rse <- create_rse_manual2(project = project_info$project, 
                           project_home = project_info$project_home, type = type, 
                           organism = project_info$organism, annotation = annotation, 
                           bfc = bfc, jxn_format = jxn_format, recount3_url = recount3_url, 
                           verbose = verbose)
  return(rse)
}