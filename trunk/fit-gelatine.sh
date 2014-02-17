FIT=./sdk/bin/fitGelFront
CUT=1800
for f in fronts/gelatine/gelatine*.dat; do
echo $f;
$FIT $f $CUT
done;

