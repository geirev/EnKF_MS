for i in *MDA-100-Adaptive-???
do
   nr=$(echo $i | cut -f1 -d'-')
   ens=$(echo $i | cut -f3 -d'-')
   trunc=$(echo $i | cut -f5 -d'-')

   name="Adaptive-1.0-${ens}-${trunc}-${nr}.dat"
   cp "${i}/rmset.dat"  ../Runs/"${name}"
done
