# Writes netcoin/colors.json: the 657 colour names of R (colors()) with their hex value.
# Run from the repository root:  Rscript tools/make_colors.R
nm <- colors()
rgbm <- col2rgb(nm)
hex <- tolower(rgb(rgbm[1, ], rgbm[2, ], rgbm[3, ], maxColorValue = 255))
names(hex) <- tolower(gsub(" ", "", nm))
hex <- hex[!duplicated(names(hex))]
body <- paste0('  "', names(hex), '": "', hex, '"', collapse = ",\n")
writeLines(c("{", body, "}"), "netcoin/colors.json")
cat("netcoin/colors.json:", length(hex), "colours\n")
