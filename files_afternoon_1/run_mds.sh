steps="nvt npt_ber_posres"
startgro="solvated_full_system.gro"

#first we need a  dummy tpr
#gmx grompp -f nvt.mdp -c $startgro -p full_system.top -o dummy.tpr -maxwarn 1
#cp full_system.top full_system_ions.top
#echo "SOL" | gmx genion -s dummy.tpr  -o full_system_ions.gro -p full_system_ions.top  -neutral  yes -conc 0.050

startgro="full_system_ions.gro"
#read -p "works ?"


for part in $steps
do
gmx grompp -f $part.mdp -c $startgro -p  full_system_ions.top -o $part.tpr -r $startgro

startgro=$part.gro
gmx mdrun -deffnm $part -v -cpi -gpu_id 1 -nt 8
done

steps="npt_ber npt" # no more position restraints here!
for part in $steps
do
gmx grompp -f $part.mdp -c $startgro -p  full_system_ions.top -o $part.tpr 
startgro=$part.gro
gmx mdrun -deffnm $part -v -cpi -gpu_id 0 -nt 8
done
