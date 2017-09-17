x = IS(f, 100, 4, 2)
z = t(apply(x, 1, function(row)which(row > 0)))

m = initMC(c(1,2))
perms = allPerm(m)
