# The same tour in R (netCoin 2.1.49), to compare with review_phase2.py
# Run from the repository root:  Rscript python_examples/review_phase2.R
suppressMessages(library(netCoin))
data(Galapagos); galapagos <- Galapagos
options(width = 160)

cat("\n== 1. dichotomize\n")
books <- data.frame(keywords = c("red;blue", "blue;green", "red", "Green;red;blue", NA, "red;red"))
print(dichotomize(books, "keywords", sep = ";"))
print(dichotomize(books, "keywords", sep = ";", nas = NULL, length = 2))
words <- data.frame(s = c("The cat and the dog", "a Dog, a cat!", "the bird"))
print(dichotomize(words, "s", sep = " ", text = "lower", stopwords = c("the", "a", "and"), nas = NULL))

cat("\n== 2. coin\n")
C <- coin(galapagos)
print(C)
print(head(diag(C), 3))
print(colnames(coin(galapagos, minimum = 10, sort = TRUE)))

cat("\n== 3. sim\n")
print(sim(C, "Jaccard")[1:3, 1:3])
print(names(sim(C, c("Jaccard", "Haberman", "confidence"))))

cat("\n== 4. edgeList\n")
print(head(edgeList(C, c("Jaccard", "Haberman"), criteria = "Jaccard", min = .7, sort = "Jaccard"), 6))
print(head(edgeList(C), 3))
print(head(expectedList(C), 3))

cat("\n== 5. coocur, propCoin, asNodes\n")
counts <- data.frame(a = c(2, 0, 1, 3, 0, 1), b = c(1, 1, 0, 2, 2, 0), c = c(0, 2, 2, 1, 1, 1))
print(coocur(counts))
print(propCoin(C, margin = 1)[1:3, 1:3])
print(head(asNodes(C, TRUE, TRUE, "es"), 3))

cat("\n== 6. allNet\n")
finches <- read.delim("python_examples/sanderson/finches.tsv")
net <- allNet(galapagos, procedures = c("Jaccard", "Haberman"), lwidth = "Jaccard", criteria = "Jaccard",
              minL = .5, nodes = finches[, c("name", "type")], color = "type", layout = "pc",
              main = "Galapagos finches (Jaccard >= 0.5)")
print(head(net$nodes, 3)); print(head(net$links, 3))
