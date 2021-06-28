## examples
## concorde and LKH2 can both give optimal orders, so just need to run one of them.

source("junli-genetic-map-functions.R")

# for Windows
concorde_path = "./bin/concorde.exe"
LKHexec = "./bin/LKH2.exe"

# for linux
concorde_path = "./bin/concorde-Linux"
LKHexec = "./bin/LKH2-Linux"

# Step 1: read your current genetic maps
m1 = read.cross(file = "rqtl-input-example.csv", format = "csvr", crosstype = "riself")

# Step2: drop markers that are too close (optional)
# for QTL mapping, high-density maps are not necessary
# it will not increase too much of accuracy of the QTL locations
# but increases the QTL analysis time a lot.
# so I suggest to thin the markers with function dropSimilarMarkers to drop markers that have a small recombination fraction

m2 = dropSimilarMarkers(m1) # drop markers with recombination fraction < 0.1 from each other
m2 = quickEst(m2) # estimate the genetic distance
m3 = tspOrder2(cross = m2, method="concorde", execPath=concorde_path) # reorder markers with tsp software concorde
# m3 = tspOrder2(cross = m2, method="lkh", execPath=LKHexec) # reorder markers with tsp software LKH2
m3 = quickEst(m3)          # estimate the genetic distance
m3 = matchOrientation(m3, m2) # make all chromosomes of m3 have the same orientation as m2

map3 = pull.map(m3, as.table = T) # pull genetic maps of m3
write.table(map3, "map3.txt",sep="\t") # write for you to compare with your current map

# Step 4: make marker orders are more consistent with the physical map (optinal)
### run multiple tspOrder2 to get the best order that is consistent
### with the physical map but without affecting the genetic map length
newmap = map3

# You can run multiple times of the loop or change the loop number (default 10)
# each run will improve the order in newmap
for (n in 1:10){
  m6 = tspOrder2(cross = m3, method="concorde", execPath=concorde_path)
  # m6 = tspOrder2(cross = m3, method="lkh", execPath=LKHexec)
  m6 = quickEst(m6)
  m6 = matchOrientation(m6, m3)
  map6 = pull.map(m6, as.table = T)
  runmaps[[(length(runmaps)+1)]] = map6
  newmap = better.order(newmap, map6)
}

# replace maps of m3 with the adjusted maps
markerlist = split(rownames(newmap), newmap$chr)
m4 = reorder.marker(m3, markerlist) # replace the maps of m3 with the new marker list
plotMap(m3, m4)
chrlen(m3)-chrlen(m4)

map4 = pull.map(m4, as.table = T) # final adjusted maps
write.table(map4, "map4.txt",sep="\t")

