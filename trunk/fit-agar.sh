FIT=./sdk/bin/fitGelFront
CUT=600
for f in fronts/agar/agar*.dat; do
echo $f;
$FIT $f $CUT
done;

