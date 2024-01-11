data(PRborder)

k <- 12

data(PRprec)
coords <- as.matrix(PRprec[sample(1:nrow(PRprec)), 1:2])

params <- c(variance = 1, kappa = 1) 

set.seed(1)
x.k <- book.rspde(coords, range = sqrt(8) / params[2], 
                  sigma = sqrt(params[1]), n = k, mesh = prmesh1,
                  return.attributes = TRUE)
