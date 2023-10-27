
case=("KS-ES")
rm -f rmse-${case}.dat tmpfile
experiments=("KS-ES-2-2" "KS-ES-3-2" "KS-ES-4-2" "KS-ES-5-2" "KS-ES-6-2" "KS-ES-7-2" "KS-ES-8-2" "KS-ES-9-2")
for ex in ${experiments[@]}
do
   echo "Exp : ${case} ${ex}"
   echo ${ex}/rmse.dat > inp.file
   a.out < inp.file >>  tmpfile
done
awk -v ln=2 '{print ln++  " "  $0 }' tmpfile > rmse-${case}.dat

case=("KS-ESX")
rm -f rmse-${case}.dat tmpfile
experiments=("KS-ES-2-2X" "KS-ES-3-2X" "KS-ES-4-2X" "KS-ES-5-2X" "KS-ES-6-2X" "KS-ES-7-2X" "KS-ES-8-2X" "KS-ES-9-2X")
for ex in ${experiments[@]}
do
   echo "Exp : ${case} ${ex}"
   echo ${ex}/rmse.dat > inp.file
   a.out < inp.file >> tmpfile
done
awk -v ln=2 '{print ln++  " "  $0 }' tmpfile > rmse-${case}.dat


case=("KS-MDA")
rm -f rmse-${case}.dat tmpfile
experiments=("KS-MDA-2-2" "KS-MDA-3-2" "KS-MDA-4-2" "KS-MDA-5-2" "KS-MDA-6-2" "KS-MDA-7-2" "KS-MDA-8-2" "KS-MDA-9-2" "KS-MDA-10-2" "KS-MDA-11-2" "KS-MDA-12-2" "KS-MDA-13-2" "KS-MDA-14-2" "KS-MDA-15-2")
for ex in ${experiments[@]}
do
   echo "Exp : ${case} ${ex}"
   echo ${ex}/rmse.dat > inp.file
   a.out < inp.file >> tmpfile
done
awk -v ln=2 '{print ln++  " "  $0 }' tmpfile > rmse-${case}.dat

case=("KS-MDAX")
rm -f rmse-${case}.dat tmpfile
experiments=("KS-MDA-2-2X" "KS-MDA-3-2X" "KS-MDA-4-2X" "KS-MDA-5-2X" "KS-MDA-6-2X" "KS-MDA-7-2X" "KS-MDA-8-2X" "KS-MDA-9-2X" "KS-MDA-10-2X" "KS-MDA-11-2X" "KS-MDA-12-2X" "KS-MDA-13-2X" "KS-MDA-14-2X" "KS-MDA-15-2X")
for ex in ${experiments[@]}
do
   echo "Exp : ${case} ${ex}"
   echo ${ex}/rmse.dat > inp.file
   a.out < inp.file >> tmpfile
done
awk -v ln=2 '{print ln++  " "  $0 }' tmpfile > rmse-${case}.dat


case=("KS-IES")
rm -f rmse-${case}.dat tmpfile
experiments=("KS-IES-2-2" "KS-IES-3-2" "KS-IES-4-2" "KS-IES-5-2" "KS-IES-6-2" "KS-IES-7-2" "KS-IES-8-2" "KS-IES-9-2" "KS-IES-10-2" "KS-IES-11-2" "KS-IES-12-2" "KS-IES-13-2" "KS-IES-14-2" "KS-IES-15-2")
for ex in ${experiments[@]}
do
   echo "Exp : ${case} ${ex}"
   echo ${ex}/rmse.dat > inp.file
   a.out < inp.file >> tmpfile
done
awk -v ln=2 '{print ln++  " "  $0 }' tmpfile > rmse-${case}.dat

case=("KS-IES4X")
rm -f rmse-${case}.dat tmpfile
experiments=("KS-IES4-2-2X" "KS-IES4-3-2X" "KS-IES4-4-2X" "KS-IES4-5-2X" "KS-IES4-6-2X" "KS-IES4-7-2X" "KS-IES4-8-2X" "KS-IES4-9-2X" "KS-IES4-10-2X" "KS-IES4-11-2X" "KS-IES4-12-2X" "KS-IES4-13-2X" "KS-IES4-14-2X" "KS-IES4-15-2X")
for ex in ${experiments[@]}
do
   echo "Exp : ${case} ${ex}"
   echo ${ex}/rmse.dat > inp.file
   a.out < inp.file >> tmpfile
done
awk -v ln=2 '{print ln++  " "  $0 }' tmpfile > rmse-${case}.dat

case=("KS-IES5X")
rm -f rmse-${case}.dat tmpfile
experiments=("KS-IES5-2-2X" "KS-IES5-3-2X" "KS-IES5-4-2X" "KS-IES5-5-2X" "KS-IES5-6-2X" "KS-IES5-7-2X" "KS-IES5-8-2X" "KS-IES5-9-2X" "KS-IES5-10-2X" "KS-IES5-11-2X" "KS-IES5-12-2X" "KS-IES5-13-2X" "KS-IES5-14-2X" "KS-IES5-15-2X")
for ex in ${experiments[@]}
do
   echo "Exp : ${case} ${ex}"
   echo ${ex}/rmse.dat > inp.file
   a.out < inp.file >> tmpfile
done
awk -v ln=2 '{print ln++  " "  $0 }' tmpfile > rmse-${case}.dat

case=("KS-IESX")
rm -f rmse-${case}.dat tmpfile
experiments=("KS-IES-2-2X" "KS-IES-3-2X" "KS-IES-4-2X" "KS-IES-5-2X" "KS-IES-6-2X" "KS-IES-7-2X" "KS-IES-8-2X" "KS-IES-9-2X" "KS-IES-10-2X" "KS-IES-11-2X" "KS-IES-12-2X" "KS-IES-13-2X" "KS-IES-14-2X" "KS-IES-15-2X")
for ex in ${experiments[@]}
do
   echo "Exp : ${case} ${ex}"
   echo ${ex}/rmse.dat > inp.file
   a.out < inp.file >> tmpfile
done
awk -v ln=2 '{print ln++  " "  $0 }' tmpfile > rmse-${case}.dat


case=("KS-MDAstep5")
rm -f rmse-${case}.dat tmpfile
experiments=("KS-MDA-5-2-1" "KS-MDA-5-2-2" "KS-MDA-5-2-3" "KS-MDA-5-2-4" "KS-MDA-5-2-5" "KS-MDA-5-2-6" "KS-MDA-5-2-7" "KS-MDA-5-2-8" "KS-MDA-5-2-9")
for ex in ${experiments[@]}
do
   echo "Exp : ${case} ${ex}"
   echo ${ex}/rmse.dat > inp.file
   a.out < inp.file >> tmpfile
done
awk -v ln=1 '{print ln++  " "  $0 }' tmpfile > rmse-${case}.dat

case=("KS-MDAstep10")
rm -f rmse-${case}.dat tmpfile
experiments=("KS-MDA-10-2-1" "KS-MDA-10-2-2" "KS-MDA-10-2-3" "KS-MDA-10-2-4" "KS-MDA-10-2-5" "KS-MDA-10-2-6" "KS-MDA-10-2-7" "KS-MDA-10-2-8" "KS-MDA-10-2-9")
for ex in ${experiments[@]}
do
   echo "Exp : ${case} ${ex}"
   echo ${ex}/rmse.dat > inp.file
   a.out < inp.file >> tmpfile
done
awk -v ln=1 '{print ln++  " "  $0 }' tmpfile > rmse-${case}.dat

case=("KS-MDAstep15")
rm -f rmse-${case}.dat tmpfile
experiments=("KS-MDA-15-2-1" "KS-MDA-15-2-2" "KS-MDA-15-2-3" "KS-MDA-15-2-4" "KS-MDA-15-2-5" "KS-MDA-15-2-6" "KS-MDA-15-2-7" "KS-MDA-15-2-8"
"KS-MDA-15-2-9" "KS-MDA-15-2-10" "KS-MDA-15-2-11")
for ex in ${experiments[@]}
do
   echo "Exp : ${case} ${ex}"
   echo ${ex}/rmse.dat > inp.file
   a.out < inp.file >> tmpfile
done
awk -v ln=1 '{print ln++  " "  $0 }' tmpfile > rmse-${case}.dat
