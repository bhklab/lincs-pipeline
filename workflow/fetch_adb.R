library(httr2)
library(jsonlite)
library(tibble)

adb <- request("https://v2annotationdb.bhklab.ca/compound/all") |>
  req_perform() |>
  resp_body_string() |>
  fromJSON(flatten = TRUE) |>
  as_tibble()

write.csv(adb,"../data/rawdata/all_adb_compounds.csv")