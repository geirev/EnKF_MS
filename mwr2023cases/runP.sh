export OMP_NUM_THREADS=14
for i in infile.KS-ES-?-2\
         infile.KS-ES-?-2X\
         infile.KS-MDA-*-2\
         infile.KS-MDA-*-2X\
         infile.KS-IES-*-2\
         infile.KS-IES-*-2X\
         infile.KS-MDA-12-5\
         infile.KS-MDA-12-5X\
         infile.PRED?\
         infile.KS-MDA-5 infile.KS-MDA-5sep infile.KS-MDA-5*-[AO]*\
         infile.KS-MDA-5-2-? infile.KS-MDA-15-2-?  infile.KS-MDA-10-2-?
for i in infile.KS-ES-?-2\
do
   echo $i
   sep=$(echo ${i:(-3)})
   exp=$(echo "${i/infile./}")
   expinfile=$(sed -n '2p' $i | sed -e 's/  */ /g' | cut -f2 -d' ')
   [[ $exp != $expinfile ]] && sed -i "s/$expinfile/$exp/g" $i

   if [ -f $i ]
   then
       cp $i infile.in
       pushd $exp
#       gnuplot ../cpdf2.gnu
       gnuplot ../rms.gnu
#       [[ -f costf.dat ]] &&  gnuplot ../costf.gnu
       popd
   fi
done
