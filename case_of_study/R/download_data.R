#' Download GBIF Occurrences
#' 
#' Downloads GBIF occurrence data for given species names.
#' It retrieves taxon keys from scientific names, downloads data on the fly,
#' and optionally processes and saves the data.
#' 
#' @param sp_names A character vector of species scientific names.
#' @param user GBIF user name for authentication.
#' @param pwd GBIF password for authentication.
#' @param email GBIF email for authentication.
#' @param start_year Starting year for occurrence data (default is NULL, it will download from the first record).
#' @param end_year Ending year for occurrence data (default is NULL, it will download till the last record).
#' @param process_data Logical indicating whether to process the downloaded data (default is TRUE).
#' @param path Path to save the downloaded data (optional).
#' @return A list of dataframes containing downloaded GBIF occurrence data.
#' @export
download_gbif_occurrences <- function(sp_names, user, pwd, email, start_year = NULL, end_year = NULL, process_data = TRUE, path = NULL) {
  message("Running download_gbif_occurrences")
  
  # Function to retrieve taxon key from scientific name
  get_taxon_key <- function(name) {
    result <- rgbif::occ_search(scientificName = name, hasCoordinate = TRUE)
    if (!is.null(result$data)) {
      return(as.character(result$data$taxonKey[1]))
    } else {
      warning(paste("No taxon key found for", name))
      return(NA)
    }
  }
  
  # Validate input parameters
  stopifnot(length(sp_names) > 0)
  stopifnot(is.character(user), is.character(pwd), is.character(email))
  stopifnot(is.null(start_year) || is.numeric(start_year))
  stopifnot(is.null(end_year) || is.numeric(end_year))
  stopifnot(is.logical(process_data))
  stopifnot(is.null(path) || is.character(path))
  
  # Transform scientific names to GBIF taxon key
  taxonKeys <- purrr::map_chr(sp_names, get_taxon_key)
  
  # Download data on the fly
  gbif_occurrences <- list()
  gbif_references <- list()
  
  for (i in seq_along(sp_names)) {
    # Obtain the taxon key for the current species
    taxon_key <- taxonKeys[i]
    
    # Initiate GBIF data and reference lists for the current species
    data_gbif <- list()
    
    # Construct predicate list
    predicates <- list(
      rgbif::pred("taxonKey", taxon_key),
      rgbif::pred("hasCoordinate", TRUE),
      rgbif::pred_gt("decimalLatitude", -90),
      rgbif::pred_gt("decimalLongitude", -180),
     rgbif::pred("occurrenceStatus", "PRESENT")
    )
    
    if (!is.null(start_year)) {
      predicates[[length(predicates)+1]] <-  rgbif::pred_gte("YEAR", start_year)
    }
      
    if (!is.null(end_year)) {
      predicates[[length(predicates)+1]] <-  rgbif::pred_lte("YEAR", end_year)
    }
    
    predicates <- do.call(rgbif::pred_and, predicates)
    
    # Download GBIF occurrence data
    request_id <- as.character(rgbif::occ_download(
      predicates,
      format = "SIMPLE_CSV",
      user = user,
      pwd = pwd,
      email = email
    ))
    
    # Wait for the download to complete
    response <- rgbif::occ_download_wait(request_id)
    
    if (response$status == "SUCCEEDED") {
      # Temporary file to store the downloaded data
      temp <- tempfile(fileext = ".csv")
      
      # Download and read the data
      download.file(response$downloadLink, temp, mode = "wb")
      data_gbif <- read.csv(unz(temp, paste0(response$key, ".csv")),
                            header = TRUE, sep = "\t", dec = ".")
      
      # Process the data if required
      if (process_data) {
        # Select relevant columns, arrange by year, and drop NAs
        data_gbif <- data_gbif %>%
          dplyr::select(decimalLongitude, decimalLatitude, year, species, countryCode, gbifID, coordinateUncertaintyInMeters) %>%
          dplyr::arrange(year) %>%
          dplyr::filter(!is.na(decimalLongitude)) %>%
          dplyr::filter(!is.na(decimalLatitude))
      }
      
      # Save or store the processed data based on the 'path' parameter
      if (!is.null(path)) {
        write.table(data_gbif,
                    file = file.path(path, paste0(sub(" ", "_", sp_names[i]), ".csv")),
                    quote = FALSE, sep = "\t", dec = ".",
                    row.names = FALSE, col.names = TRUE)
      } else {
        gbif_occurrences[[sp_names[i]]] <- data_gbif
      }
      
      # Store the GBIF reference link
      gbif_references[[sp_names[i]]] <- rgbif::gbif_citation(request_id)$download
      
      # Delete the temporary file
      unlink(temp)
      
    } else {
      # Issue a warning if the download fails
      warning(paste("Download for", sp_names[i], "failed:", response$status))
    }
  }
  
  # Display the GBIF reference links
  for (name in names(gbif_references)) {
    message(paste0(name, ": ", gbif_references[[name]]))
  }
  
  # Return the processed data
  return(gbif_occurrences)
}
