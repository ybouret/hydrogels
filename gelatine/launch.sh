../sdk/bin/rdiff1d -f "Cgel=0;"    ../src/rdiff1d/gelatine.lua ../src/rdiff1d/setup.lua 
for p in 8 \
7.9 7.8 7.7 7.6 7.5 7.4 7.3 7.2 7.1 7.0 \
6.9 6.8 6.7 6.6 6.5 6.4 6.3 6.2 6.1 6.0 \
5.9 5.8 5.7 5.6 5.5 5.4 5.3 5.2 5.1 5.0 \
4.9 4.8 4.7 4.6 4.5 4.4 4.3 4.2 4.1 4.0 \
3.9 3.8 3.7 3.6 3.5 3.4 3.3 3.2 3.1 3.0 \
2.9 2.8 2.7 2.6 2.5 2.4 2.3 2.2 2.1 2.0 \
1.9 1.8 1.7 1.6 1.5 1.4 1.3 1.2 1.1 1.0;
do
../sdk/bin/rdiff1d -f "Cgel=10^(-$p);" ../src/rdiff1d/gelatine.lua ../src/rdiff1d/setup.lua
done