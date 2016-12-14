#!/usr/bin/env Rscript
library(pipeR)
library(tidyverse)
library(rgl)

library(tumorr)

#########1#########2#########3#########4#########5#########6#########7#########
## transformation

if (rgl.cur()) {rgl.close()}
rgl::open3d(windowRect=c(0, 0, 600, 600))
rgl::clear3d()
rgl::view3d(70, 5, 60)
axes3d()
.r = 2
seq(-.r, .r) %>>%
    (expand.grid(x=., y=., z=.)) %>>%
    tumorr:::trans_coord_hcc() %>>%
    # tumorr:::trans_coord_fcc() %>>%
    dplyr::mutate(r= sqrt(x*x+y*y+z*z)) %>>%
    with(spheres3d(x, y, z, color='#009999', radius=0.51, alpha=0.6))
title3d('', '', 'x', 'y', 'z')

#########1#########2#########3#########4#########5#########6#########7#########
## minimum

.hex_xy = read_csv('x,y,z
0,0,0
1,0,0
0,1,0
1,0,-1')

if (rgl.cur()) {rgl.close()}
rgl::open3d(windowRect=c(0, 0, 600, 600))
rgl::clear3d()
rgl::view3d(20, 10, 60)
axes3d()
.hex_xy %>>%
    trans_coord_hex() %>>%
    with(spheres3d(x, y, z, color='#009999', radius=0.51, alpha=0.6))
title3d('', '', 'x', 'y', 'z')

#########1#########2#########3#########4#########5#########6#########7#########
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
axes3d()
.hex_xy %>>%
    bind_rows(.hex_xy %>>% dplyr::filter(x > 0 | (x==0 & y==0)) %>>% mutate(z=-1)) %>>%
    bind_rows(.hex_xy %>>% dplyr::filter(x < 0 | (x==0 & y==0)) %>>% mutate(z=1)) %>>% (?.) %>>%
    tumorr:::trans_coord_fcc() %>>%
    with(spheres3d(x, y, z, color='#009999', radius=0.51, alpha=0.6))
title3d('', '', 'x', 'y', 'z')


if (rgl.cur()) {rgl.close()}
rgl::open3d(windowRect=c(0, 0, 600, 600))
rgl::clear3d()
rgl::view3d(40, 20, 60)
axes3d()
.hex_xy %>>%
    bind_rows(.hex_xy %>>% mutate(z=-1) %>>% dplyr::filter(x < 0 | (x==0 & y==0))) %>>%
    bind_rows(.hex_xy %>>% mutate(z=1) %>>% dplyr::filter(x < 0 | (x==0 & y==0))) %>>% (?.) %>>%
    bind_rows(.hex_xy %>>% mutate(z=5)) %>>%
    bind_rows(.hex_xy %>>% mutate(z=4) %>>% dplyr::filter(x > 0 | (x==0 & y==0))) %>>%
    bind_rows(.hex_xy %>>% mutate(z=6) %>>% dplyr::filter(x > 0 | (x==0 & y==0))) %>>%
    tumorr:::trans_coord_hcc() %>>%
    with(spheres3d(x, y, z, color='#009999', radius=0.51, alpha=0.6))
title3d('', '', 'x', 'y', 'z')
