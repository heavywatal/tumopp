library(pipeR)
library(tidyverse)
library(wtl)
library(tumorr)
library(igraph)

refresh('tumorr')
(.result = tumopp(str_split('-D3 -Chex -Pifany', ' ')[[1]]))
(.population = .result$population[[1]])
(.extant = .population %>>% filter_extant())

.center = dplyr::sample_n(.extant %>>% dplyr::filter(z == 0), 1) %>>% (?.)
.specimen = sample_bulk(.extant, .center)

.center2 = dplyr::sample_n(.extant %>>% dplyr::filter(z == 0), 1) %>>% (?.)
.specimen2 = sample_bulk(.extant, .center2)

gglattice2d(.extant %>>% dplyr::filter(-1 <= z, z < 1), 'z', , 0.5)+
geom_point(data=.specimen %>>% dplyr::filter(-1 <= z, z < 1), colour='red')+
geom_point(data=.specimen2 %>>% dplyr::filter(-1 <= z, z < 1), colour='blue')

.el = .population %>>%
#    filter_connected(.specimen$genealogy) %>>%
    make_edgelist()

.g = .el %>>% igraph::graph_from_data_frame()

.el %>>% layout_genealogy() %>>% ggplot_genealogy() %>>% grid::grid.draw()

# T_S
(.ts1 = mean_branch_length(.g, as.character(.specimen$id)))
(.ts2 = mean_branch_length(.g, as.character(.specimen2$id)))
(.ts = mean(.ts1, .ts2))

# T_B
(.tb = mean_branch_length(.g, as.character(.specimen$id), as.character(.specimen2$id)))

# T_D
(.td = .tb - .ts)

# T_T
(.tt = (.tb + .ts) / 2)

mean_branch_length(.g, c(as.character(.specimen$id), as.character(.specimen2$id)))

fst_HSM(.ts, .tb)

fst_HBK(.ts, .tb)
