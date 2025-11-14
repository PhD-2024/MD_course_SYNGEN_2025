steps="steep nvt npt_ber npt"
startgro="solvated_full_system.gro"
for part in $steps
do
gmx grompp -f $part.mdp -c $startgro -p  full_system.top -o $part.tpr
read -p works?
startgro=$part.gro
gmx mdrun -deffnm $part -v 
done

