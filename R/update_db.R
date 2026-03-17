
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
