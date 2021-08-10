## examples for outcrossing crops
## data from onemap package: https://github.com/augusto-garcia/fullsibQTL
## code is borrowed from: 
# https://cran.r-project.org/web/packages/onemap/vignettes/Outcrossing_Populations.html
## concorde has some problem for some LGs, so I used LKH2 below.

source("junli-genetic-map-functions.R")

# for Windows
concorde_path = "./bin/concorde.exe"
LKHexec = "./bin/LKH2.exe"

# for linux
concorde_path = "./bin/concorde-Linux"
LKHexec = "./bin/LKH2-Linux"


# Step 1: read your current genetic maps
library(onemap)
onemap_example_out = read_onemap(inputfile = "example_QTLfullsib.raw")
plot(onemap_example_out)
plot_by_segreg_type(onemap_example_out)

# Step 2: remove redundant markers
bins <- find_bins(onemap_example_out, exact = FALSE)
bins

# create a new onemap object without redundant markers
# no redundant markers for this dataset
bins_example <- create_data_bins(onemap_example_out, bins)
bins_example

## Step 3: Testing segregation pattern
# test_segregation_of_a_marker(bins_example, 4)
segreg_test <- test_segregation(bins_example)
# print(segreg_test)
select_segreg(segreg_test, distorted = TRUE) #to show the markers names with segregation 
# select_segreg(segreg_test, distorted = FALSE) #to show the markers names without segregation distortion

dist <- select_segreg(segreg_test, distorted = TRUE, numbers = TRUE) #to show the markers numbers with segregation distortion
dist
#> [1] 17 19 25

no_dist <- select_segreg(segreg_test, distorted = FALSE, numbers = TRUE) #to show the markers numbers without segregation distortion
# no_dist

plot(segreg_test)


## Step 4: Estimating two-point recombination fractions
LOD_sug <- suggest_lod(bins_example)
LOD_sug
twopts <- rf_2pts(bins_example, LOD = LOD_sug, max.rf = 0.5)
twopts
print(twopts, c("M1", "M3"))

## Step 5: Assigning markers to linkage groups
mark_no_dist <- make_seq(twopts, c(no_dist)) # without distorted markers
LGs <- group(mark_no_dist, LOD=4, max.rf=0.35) # change LOD criteria to split them into reasonable chromosomes
LGs
str(LGs)

## Step 6: make maps
set_map_fun(type = "kosambi")

## one loop for all LGs
ordered.LGs = list()
tol = 1e-4
for (x in 1:LGs$n.groups) {
  cat("x is", x, "\n")
  LGx <- make_seq(LGs, x)
  if (length(LGx$seq.num) < 3) next # pass if only has 1 or 2 markers
  # use TSP software
  # I found some D type markers has rf=1e-50 and LOD=0
  # I suppose they should be non-linked, so I set min.LOD=0.1
  # the tutorial set min.LOD=0, which will treat the two markers as totally linked.
  # Correct me if you know more about outcrossing species
  rfx = onemap:::get_mat_rf_out(LGx, LOD=FALSE, max.rf=0.5, min.LOD=0.1)
  rfx[is.na(rfx)] = 0.5
  diag(rfx) = NA
  # print(rfx)
  orderx = tspOrder2.rf(rfx, method="lkh", execPath=LKHexec)
  # orderx = tspOrder2.rf(rfx, method="concorde", execPath=concorde_path)
  print(orderx)
  LGx_lkh = map(make_seq(get(LGx$twopt), LGx$seq.num[onemap:::avoid_reverse(orderx)], twopt = LGx$twopt), tol = tol)
  ordered.LGs[[x]] = LGx_lkh
}
# print LG1 map and parent phase
ordered.LGs[[1]]
# print all
ordered.LGs

# draw map
draw_map(ordered.LGs, names = TRUE, grid = TRUE, cex.mrk = 1)

# to compare with the onemap function results
# warning: if you have hundreds of markers, it might take a long time.
x = 1 # change here to see different LGs
ordered.LGs[[x]]

(order.ser <- seriation(make_seq(LGs, x)))
(order.rcd <- rcd(make_seq(LGs, x)))
(order.re <- record(make_seq(LGs, x)))
(order.ug <- ug(make_seq(LGs, x)))
ordered.LGs[[x]]
# print the rf matrix
onemap:::get_mat_rf_out(make_seq(LGs, x), LOD=FALSE, max.rf=0.5, min.LOD=0)
# check a few markers
print(twopts, c("M9", "M16"))
print(twopts, c("M9", "M21"))
print(twopts, c("M23", "M21"))
print(twopts, c("M35", "M45"))
marker_type(LG2)
