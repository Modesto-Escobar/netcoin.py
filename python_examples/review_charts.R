# The same calls in R (netCoin 2.1.49), to compare with review_charts.py
# Run from the repository root:  Rscript python_examples/review_charts.R
suppressMessages(library(netCoin))
options(width = 160)
data(ess, Galapagos, sociologists, works)

cat("== barCoin\n")
# R needs the weight column among the variables (python does not)
bar <- barCoin(ess, c("Gender", "Age", "Social participation", "Political participation", "cweight"),
               dichotomies = c("Social participation", "Political participation"), valueDicho = "Yes",
               weight = "cweight", expected = TRUE, confidence = TRUE, significance = TRUE,
               select = "Gender:Male")
print(head(bar$nodes, 4)); print(head(bar$links, 3)); str(bar$options)

cat("\n== timeCoin\n")
soc <- as.data.frame(sociologists); wrk <- as.data.frame(works)
timeline <- timeCoin(soc, start = "birth", end = "death", group = "school", text = "label",
                     main = "Sociologists", events = wrk, eventNames = "works", eventPeriod = "name",
                     eventTime = "date", eventColor = "label")
print(head(timeline$periods[, c("name", "birth", "death", "school")], 3))
summary(timeline)

cat("\n== incTime / mobileEdges\n")
print(incTime(soc, "name", "birth", "death")[1:3, 1:4])
print(head(mobileEdges(soc, "name", "birth", 10), 3))

cat("\n== pieCoin\n")
pie <- pieCoin(coin(Galapagos), main = "Galapagos finches")
print(dim(pie$data)); print(pie$data[1, 2, ])
