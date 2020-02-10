#!/bin/bash

rm -f plotDat.pg
touch plotDat.pg
cat << EOF >> plotDat.pg
#!/usr/bin/gnuplot
reset
set terminal png
set xrange [-1:1]
set yrange [-1:1]
set key off
unset colorbox
set size ratio -1
set palette maxcolors 2
set palette defined ( 1 '#FF0000',\
                      2 '#FFFFFF')
plot "matlab.dat" using ((\$1-25)/25):((\$2-25)/25):3 matrix with image
EOF

chmod +x plotDat.pg
./plotDat.pg > Iterations.png
rm -f plotDat.pg
