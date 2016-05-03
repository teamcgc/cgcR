#!/usr/local/bin/Rscript
'Usage: get_http_file.R URL FILE

Options:
'  -> doc

library(docopt)
library(httr)
opts <- docopt(doc)
print(opts)

ofname = NULL

if(is.null(opts$FILE)) {
    ofname=basename(opts$URL)
} else {
    ofname=opts$FILE
}

GET(opts$URL,write_disk(ofname))
