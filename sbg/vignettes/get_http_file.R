#!/usr/local/bin/Rscript
'Usage: get_http_file.R URL FILE

Options:
'  -> doc

library(docopt)
library(httr)
opts <- docopt(doc)
print(opts)

dir.create('output')

ofname = NULL

if(is.null(opts$FILE)) {
    ofname=basename(opts$URL)
} else {
    ofname=opts$FILE
}

ofname=file.path('output',ofname)

GET(opts$URL,write_disk(ofname))
