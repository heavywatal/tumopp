#!/usr/bin/env Rscript
library(rgl)

library(tumorr)


## transformation

if (rgl.cur()) {rgl.close()}
rgl::open3d(windowRect=c(0, 0, 600, 600))
rgl::clear3d()
rgl::view3d(70, 5, 60)
clear3d()
axes3d()
.r = 2
seq(-.r, .r) %>>%
    (expand.grid(x=., y=., z=.)) %>>%
#    trans_coord_hcc %>>%
    trans_coord_fcc %>>%
    mutate(r= sqrt(x*x+y*y+z*z)) %>>%
#    dplyr::filter(r <= 3) %>>% (?.) %>>%
with(spheres3d(x, y, z, color='#009999', radius=0.51, alpha=0.6))
title3d('', '', 'x', 'y', 'z')


## cells, volume, radius

.r = 60
.cells = seq(-.r, .r) %>>%
    (expand.grid(grid='regular', x=., y=., z=.)) %>>%
    tbl_df %>>%
    bind_rows(trans_coord_fcc(.) %>>% mutate(grid='hex')) %>>%
    mutate(r= sqrt(x*x+y*y+z*z)) %>>%
    dplyr::filter(r < 0.75 * .r) %>>%
    group_by(grid) %>>%
    mutate(n = order(order(r))) %>>%
    (?.)

.curve = data_frame(grid='regular', r=seq_len(100) / 2) %>>%
    mutate(n=sphere_volume(r)) %>>%
    bind_rows(mutate(., grid='hex', n= sqrt(2) * n)) %>>%
    (?.)

.p = .cells %>>%
    ggplot(aes(n, r, colour=grid))+
    geom_point()+
    geom_path(data=.curve)
.p
ggsave('radius-nodes.png', .p, width=7, height=7)

.cells %>>% group_by(grid) %>>% tally
summary(.cells)

if (rgl.cur()) {rgl.close()}
rgl::open3d(windowRect=c(0, 0, 600, 600))
rgl::clear3d()
rgl::view3d(70, 5, 60)
clear3d()
axes3d()
.cells %>>%
    dplyr::filter(r > 47) %>>% (?.) %>>%
with(spheres3d(x, y, z, color='#009999', radius=0.51, alpha=0.6))
title3d('', '', 'x', 'y', 'z')


## minimum

.hex_xy = read_csv('x,y,z
0,0,0
1,0,0
0,1,0
1,0,-1
')

if (rgl.cur()) {rgl.close()}
rgl::open3d(windowRect=c(0, 0, 600, 600))
rgl::clear3d()
rgl::view3d(20, 10, 60)
clear3d()
axes3d()
.hex_xy %>>%
    trans_coord_fcc %>>%
with(spheres3d(x, y, z, color='#009999', radius=0.51, alpha=0.6))
title3d('', '', 'x', 'y', 'z')

## neighbors

.hex_xy = read_csv('x,y,z
0,0,0
0,1,0
0,-1,0
-1,0,0
-1,1,0
1,0,0
1,-1,0')

if (rgl.cur()) {rgl.close()}
rgl::open3d(windowRect=c(0, 0, 600, 600))
rgl::clear3d()
rgl::view3d(20, 10, 60)
clear3d()
axes3d()
.hex_xy %>>%
    bind_rows(.hex_xy %>>% dplyr::filter(x > 0 | (x==0 & y==0)) %>>% mutate(z=-1)) %>>%
    bind_rows(.hex_xy %>>% dplyr::filter(x < 0 | (x==0 & y==0)) %>>% mutate(z=1)) %>>% (?.) %>>%
    trans_coord_fcc %>>%
with(spheres3d(x, y, z, color='#009999', radius=0.51, alpha=0.6))
title3d('', '', 'x', 'y', 'z')


if (rgl.cur()) {rgl.close()}
rgl::open3d(windowRect=c(0, 0, 600, 600))
rgl::clear3d()
rgl::view3d(40, 20, 60)
clear3d()
axes3d()
.hex_xy %>>%
    bind_rows(.hex_xy %>>% mutate(z=-1) %>>% dplyr::filter(x < 0 | (x==0 & y==0))) %>>%
    bind_rows(.hex_xy %>>% mutate(z=1) %>>% dplyr::filter(x < 0 | (x==0 & y==0))) %>>% (?.) %>>%
    bind_rows(.hex_xy %>>% mutate(z=5)) %>>%
    bind_rows(.hex_xy %>>% mutate(z=4) %>>% dplyr::filter(x > 0 | (x==0 & y==0))) %>>%
    bind_rows(.hex_xy %>>% mutate(z=6) %>>% dplyr::filter(x > 0 | (x==0 & y==0))) %>>%
    trans_coord_hcc %>>%
with(spheres3d(x, y, z, color='#009999', radius=0.51, alpha=0.6))
title3d('', '', 'x', 'y', 'z')
