## examples
# library(ASMap)
source("junli-genetic-map-functions.R")
# for windows
concorde_path = "./bin/concorde.exe"
LKHexec = "./bin/LKH2.exe"

# for linux
concorde_path = "./bin/concorde-Linux"
LKHexec = "./bin/LKH2-Linux"

m1 = read.cross(file = "rqtl-input-example.csv", format = "csvr", crosstype = "riself")

# drop similar markers
m2 = dropSimilarMarkers(m1)
m2 = quickEst(m2)
system.time({ m3 = tspOrder2(cross = m2, method="concorde", execPath=concorde_path) })
m3 = quickEst(m3)
m3 = matchOrientation(m3, m2)

system.time({ m4 = tspOrder2(cross = m2, method="lkh", execPath=LKHexec) })
m4 = quickEst(m4)
m4 = matchOrientation(m4, m3)

map3 = pull.map(m3, as.table = T)
write.table(map3, "map3.txt",sep="\t")
map4 = pull.map(m4, as.table = T)
write.table(map4, "map4.txt",sep="\t")


### run multiple tspOrder2 to get the best order that is consistent
### with the physical map but without affecting the genetic map length
newmap = map3

# You can run multiple times of the loop
# or change to more loops
runmaps = list()
for (n in 1:10){
  # m6 = tspOrder2(cross = m3, method="concorde", execPath=concorde_path)
  m6 = tspOrder2(cross = m3, method="lkh", execPath=LKHexec)
  m6 = quickEst(m6)
  m6 = matchOrientation(m6, m3)
  map6 = pull.map(m6, as.table = T)
  runmaps[[(length(runmaps)+1)]] = map6
  newmap = better.order(newmap, map6)
}

rowdiffs = lapply(runmaps, function(x){
  which(rownames(map3)!=rownames(x))
})

pairdif = c()
for (i in 1:length(rowdiffs)){
  for (j in 1:length(rowdiffs)){
    pairdif = c(pairdif, length(intersect(rowdiffs[[i]],rowdiffs[[j]]))/length(union(rowdiffs[[i]],rowdiffs[[j]])))
  }
}

matrix(pairdif,nrow=length(rowdiffs))

length(unique(unlist(rowdiffs)))

markerlist = split(rownames(newmap), newmap$chr)

m7 = reorder.marker(m3, markerlist)
plotMap(m3, m7)
chrlen(m3)-chrlen(m7)

map7 = pull.map(m7, as.table = T)
write.table(map7, "map7.txt",sep="\t")

