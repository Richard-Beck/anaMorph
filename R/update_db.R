
# Read database credentials from a simple KEY=value file.
# Input: `filepath` path to a text file containing HOST, DBNAME, USER, PASSWORD.
# Output: named character vector of credential values.
load_db_vars <- function(filepath) {
  if (!file.exists(filepath)) {
    stop(sprintf("Error: File '%s' does not exist.", filepath))
  }
  
  lines <- tryCatch(readLines(filepath, warn = FALSE), 
                    error = function(e) stop("Error: Failed to read file '", filepath, "'. ", e$message))
  
  if (length(lines) == 0) {
    stop(sprintf("Error: File '%s' is empty.", filepath))
  }
  
  if (any(!grepl("^[A-Za-z_][A-Za-z0-9_]*=.+$", lines))) {
    stop("Error: All lines must be in the format KEY=value, with valid variable names.")
  }
  
  kv_pairs <- strsplit(lines, "=", fixed = TRUE)
  keys <- vapply(kv_pairs, `[`, "", 1)
  values <- vapply(kv_pairs, `[`, "", 2)
  vars <- setNames(values, keys)
  
  required_keys <- c("HOST", "DBNAME", "USER", "PASSWORD")
  missing_keys <- setdiff(required_keys, names(vars))
  if (length(missing_keys) > 0) {
    stop("Error: Missing required keys: ", paste(missing_keys, collapse = ", "))
  }
  
  return(vars)
}

update_db <- function(){
  require(DBI)
  require(RMariaDB)
  require(data.table)
  
  dbvars <- load_db_vars("db_creds.txt")
  db <- dbConnect(
    MariaDB(),
    host = dbvars["HOST"],
    user = dbvars["USER"],
    password = dbvars["PASSWORD"],
    dbname = dbvars["DBNAME"]
  )
  
  passaging_raw <- dbGetQuery(db, "SELECT * FROM Passaging")
  media_raw <- dbGetQuery(db, "SELECT * FROM Media")
  dbDisconnect(db)
  
  fwrite(passaging_raw,"core_data/passaging.csv",sep=",")
  fwrite(media_raw,"core_data/media.csv",sep=",")
  
}
