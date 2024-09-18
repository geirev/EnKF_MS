#!/bin/bash

rundir="Runs"

for LOC in 1 2
do
   [[ ${LOC} -eq 1 ]] && XLOC="Distance"
   [[ ${LOC} -eq 2 ]] && XLOC="Adaptive"

   DAMPLIST=($(ls Runs/${XLOC}* | cut -f2 -d"-" | sort -u))
   ENSLIST=($(ls Runs/${XLOC}*  | cut -f3 -d"-" | sort -u))
   DISTLIST=($(ls Runs/${XLOC}* | cut -f4 -d"-" | sort -u))

   for n in  "${ENSLIST[@]}"
   do

      for DAMP in  "${DAMPLIST[@]}"
      do

         EX="${XLOC}"-"${DAMP}"-"${n}"
         for DIST in  "${DISTLIST[@]}"
         do

            EXP="${XLOC}"-"${DAMP}"-"${n}"-"${DIST}"
            echo Processing case "${EXP}"
            rm -f "${EXP}".dat

            num=$(find ./"${rundir}" -name "${EXP}-*" -print | wc -l)
            echo nrfiles "${num}"

            if [ "${num}" -gt 0 ]
            then
               cat Runs/"${EXP}"-* | awk -v n="${num}" '{for(i=1;i<=4;i++)$i=(a[i]+=$i/n)}END{print}'  > tot_"${EXP}".dat

               echo "${DIST}" "${n}"  > ttt
               paste -d'  ' ttt tot_"${EXP}".dat  > rms_"${EXP}".dat
               rm tot_"${EXP}".dat
            fi

         done
         nrf=$(find . -name "rms_${EX}*.dat" -print | wc -l)
         if [ "${nrf}" -gt 0 ]
         then
            echo "VARIABLES = \"Trunc\" \"nrens\" \"rmseA\" \"rmseO\" \"rmssA\" \"rmssO\"" > rms_all_"${EX}".dat
            echo "ZONE  T=\"${EX}\" F=Point, I= ${nrf}  , J=   1, K=1" >> rms_all_"${EX}".dat
            cat rms_"${EX}"*.dat | sort -n >> rms_all_"${EX}".dat
            rm rms_"${EX}"*.dat
         fi

      done
   done
done
