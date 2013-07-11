#!/bin/bash
#./MolproToMcci.sh filename 

input=${1}



if [  ! -e $input ]; then
echo $input not found
exit
fi


grep 'Frozen orbitals:' $input | sed -n 's/[^(]*(//p' | sed s/')'//g > frozmolpro.in


grep 'Active orbitals:' $input |sed -n 's/[^(]*(//p' | sed s/')'//g > activemolpro.in

frozorbs=( $( cat frozmolpro.in ) )

activeorbs=( $( cat activemolpro.in ) )



sizeactive=$(( ${#activeorbs[@]} ))



cumulativeorbs[0]=0

for i in $(eval echo {1..$(($sizeactive-1))})
do

j=$(($i -1))
cumulativeorbs[$i]=$((${activeorbs[$j]}+${cumulativeorbs[$j]}))
done









spin=`grep 'Spin quantum number:' $input | sed s/'Spin quantum number:'//g`

if grep -q 'Final alpha' $input ; then



k=$(grep 'Final alpha' $input | sed s/'Final alpha occupancy:'//g)
cat >occmolpro.in << EOF
$k 
EOF


alphaorbs=( $( cat occmolpro.in ) )

nsym=$(( ${#alphaorbs[@]} ))


# loop through and remove frozen 
for i in $(eval echo {0..$(($nsym-1))})
do
alphaorbs[$i]=$(( ${alphaorbs[$i]}-${frozorbs[$i]} ))
done





# convert molpro number of orbitals of each symmetry to orbital numbers
eu=0

for i in $(eval echo {0..$(($nsym-1))})
do

if [ "${alphaorbs[$i]}" -eq "0" ]
  then
  continue   
  fi

for j in $(eval echo {0..$((${alphaorbs[$i]}-1))})
do

ku[$eu]=$(($j+1+${cumulativeorbs[$i]}))

eu=$(($eu +1))
done
done






k=$(grep 'Final beta' $input | sed s/'Final beta  occupancy:'//g)
cat >occmolpro.in << EOF
$k 
EOF

betaorbs=( $( cat occmolpro.in ) )

nsym=$(( ${#betaorbs[@]} ))





# loop through and remove frozen 
for i in $(eval echo {0..$(($nsym-1))})
do
betaorbs[$i]=$(( ${betaorbs[$i]}-${frozorbs[$i]} ))
done



ed=0

for i in $(eval echo {0..$(($nsym-1))})
do

if [ "${betaorbs[$i]}" -eq "0" ]
  then
	continue   
  fi

for j in $(eval echo {0..$((${betaorbs[$i]}-1))})
do

kd[$ed]=$(($j+1+${cumulativeorbs[$i]}))

ed=$(($ed +1))
done
done




else


k=$(grep ' Final occupancy:' $input | sed s/' Final occupancy:'//g)

cat >occmolpro.in << EOF
$k 
EOF

alphaorbs=( $( cat occmolpro.in ) )

nsym=$(( ${#alphaorbs[@]} ))





# loop through and remove frozen 
for i in $(eval echo {0..$(($nsym-1))})
do
alphaorbs[$i]=$(( ${alphaorbs[$i]}-${frozorbs[$i]} ))
done





# convert molpro number of orbitals of each symmetry to orbital numbers
eu=0

for i in $(eval echo {0..$(($nsym-1))})
do

if [ "${alphaorbs[$i]}" -eq "0" ]
  then
	continue   
  fi

for j in $(eval echo {0..$((${alphaorbs[$i]}-1))})
do

ku[$eu]=$(($j+1+${cumulativeorbs[$i]}))

eu=$(($eu +1))
done
done








kd="${ku[@]}"
ed=$eu

fi

# Add commas to orbital lists


ku=`echo ${ku[@]} | sed 's/ /,/g'`
kd=`echo ${kd[@]} | sed 's/ /,/g'`

cat >mcci.in <<EOF
! mcci.in file for version four of mcci
restart         = .false.                 ! .TRUE. sets inflg = 1, .FALSE sets inflg = 0
n_up            = $eu                     ! number of up electrons
mo_up           = $ku! occupied spin up orbital labels
n_dn            = $ed
mo_dn           = $kd ! occupied spin up orbital labels
s               = $spin                     ! spin ( units of h_bar ) (principle case m_s = s)
ieig            = 1                       ! ith eigenvalue
frac            = 1.0
maxtry          = 3000                ! how many diagonalisations?
cmin            = 5.0d-4
frozen_doubly   =  
npfull          = 10
i_want_conv     = .true.
npfull_conv     = .true.
conv_thresh_e   = 0.001
conv_thresh_l   = 100.0
conv_average    = 3            ! 3 fullprune steps are being averaged
conv_history    = 3            ! 3 averaged Dvalues are being tracked
SCF_integral_filename= "FCIDUMP"
!parameters
maxc = 200000
kmax = 35
maxh = 400000000
EOF


rm occmolpro.in
rm frozmolpro.in
rm activemolpro.in
