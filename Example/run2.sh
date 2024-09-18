#!/bin/bash
export OMP_NUM_THREADS=8

LOC=1
ENSLIST=(100)
DAMPLIST=(1.0 2.0 4.0 8.0)
nrsamp=50

[[ ${LOC} -eq 1 ]] && XLOC="Distance"
[[ ${LOC} -eq 1 ]] && DISTLIST=(25 50 55 60 65 70 75 100 125 150)
[[ ${LOC} -eq 2 ]] && XLOC="Adaptive"
[[ ${LOC} -eq 2 ]] && DISTLIST=(2.2 2.4 2.6 2.8 3.0 3.2 3.4)

iter=0
while [ ${iter} -lt "${nrsamp}" ]
do
   (( iter++ ))

   for n in  "${ENSLIST[@]}"
   do

      for DAMP in  "${DAMPLIST[@]}"
      do

         for DIST in  "${DISTLIST[@]}"
         do

            attempt=0
            while [ ${attempt} -lt 5 ]
            do
               (( attempt++ ))
               [[ -f seed.dat ]] && rm seed.dat

               EXPERIMENT=${XLOC}-${DAMP}-${n}-${DIST}-${iter}
               INFILE=infile.${EXPERIMENT}

               sed -e "s/ENSSIZE/${n}/"  -e "s/XLOC/${LOC}/"  -e "s/XDIST/${DIST}.0/" -e "s/XDAMP/${DAMP}/"  -e "s/XTRUNC/${DIST}/"  < Xinfile.BASE > "${INFILE}"
               expinfile=$(sed -n '2p' "${INFILE}" | sed -e 's/  */ /g' | cut -f2 -d' ')
               [[ ${EXPERIMENT} != "${expinfile}" ]] && sed -i "s/${expinfile}/${EXPERIMENT}/g" infile."${EXPERIMENT}"

               if [ -f "${INFILE}" ]
               then
                  echo "${EXPERIMENT}"
                  cp "${INFILE}" infile.in
                  nice -19 multiscale   >   "${EXPERIMENT}".log
                  rm -f infile.in
                  RESULT=$?
                  if [ "${RESULT}" == 0 ]; then
                     echo succeded
                     rm "${EXPERIMENT}".log
                     pushd "${EXPERIMENT}" || exit
                     rm -f gnu_*.dat obs*.dat infile.in rmse.dat
                     mv rmset.dat ../Runs/"${EXPERIMENT}".dat
                     popd || exit
                     rm -f "${INFILE}"
                     rmdir "${EXPERIMENT}"
                     break
                  else
                    echo failed "${attempt}"
                  fi
               fi

            done # Attempt
         done    # Dist
      done       # Damp
   done          # enssize
done             # number of repetative runs
