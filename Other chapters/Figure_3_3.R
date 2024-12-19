
set.seed(123)

qtle <- mockQTLE()

assays(qtle)
qtle <- callSignificance(qtle)

# There is now an assay called 'significant'
assays(qtle)

pvalues(qtle)[9, 1:10]
# looking for a situation where there is a state < 0.05, and some states that are 
# bigger than that but smaller than 0.1

# Use feature-wise FDR correction -------------------------------------------
qtle_1 <- callSignificance(qtle, thresh=0.05, mode="simple")
qtle_2 <- callSignificance(qtle, thresh=0.05, secondThresh =0.1, mode="simple")

round(pvalues(qtle)[9, 1:10], 3)
assay(qtle_1, "significant")[9, 1:10]
assay(qtle_2, "significant")[9, 1:10]

