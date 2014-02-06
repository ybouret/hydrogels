for i in Cgel*.dat; do
echo $i;
../sdk/bin/fitFront $i 0
done

