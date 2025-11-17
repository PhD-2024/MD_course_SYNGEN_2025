# MD_course

The following is an overview for the practical side of the course.
You will find both information, files and instuctions here.

## Linux/Bash introduction (day1)

When working on a cluster there may not be a GUI. Also on your own
machine it may be a lot more efficient to use a terminal.

You can do everything you can do via a GUI by typing commands. 
Here we will do so using `bash` the currently still most used shell.

Even better you can automate tedious and repetitive processes, which
would require a lot of time to do by hand.

Here is a quick overview over important bash commands

| Command | Description |
|---|---|
| `.` | Current directory |
| `~` | Home directory |
| `..` | Parent directory (one level up) |
| `;` | end of line replacement if you want multiple commands in 1 line | 
| `cd destination` | Change directory (e.g. `cd /path/to/dir`) - not writing a destination goes to `~` |
| `mkdir -p` | Create a directory and any necessary parent directories |
| `mv` | Move or rename files or directories (`mv source target`) |
| `cp` | Copy files or directories (`cp source target`) |
| `cat` | concatenates files (`cat file1 file2 file3`) |
| `touch`| creates an empty file (or updates times) | 
| `scp` | Secure copy files between hosts (`scp file user@host:/path`) |
| `pwd` | Print working directory |
| `sed "s\|to_substitute\|replaced_by\|g" ` | Replace all occurrences of `to_substitute` with `replaced_by` using `sed` |
| `grep string files`| searches for string in files| 
| `man command`| opens the manual page for the command - RTFM (Read the very fine manual... |
| `ssh user@machine` | connects to machine via ssh for your username "user" on the remote machine | 
| `echo` | prints to stout: e.g. `echo foo bar`| 
| `chmod` | Change file permissions (mode). Examples: `chmod 755 file` or `chmod +x script.sh` |
| `chown` | Change file owner and group. Example: `chown user:group file` |


1) (Bash) Scripts are essentially nothing essentially nothing more storing commands in a text file so you can reuse them (also great to have reproducability - no arbitray behavior due to typos-). You usually start them with the so called shebang `#!/bin/bash` which tells the computer with which interpreter to execute it - here bash. `#!/bin/python3`would be a python script for example. You can of course completely forgo it and manually call the `bash`command in the shell.

2) Wildcards/placeholders `*` any chars, `?` one char.
3) Redirections `>` redirects stout to a new file.
`>>` appends to file.
`&>` redirects stout and sterr to a file.

3) With loops you can do a lot some examples

For loop over iterations:

```
for i in {1..5}; do
  echo "Iteration $i"
done
```


```
while true; do
    sleep 5
    # do something every 5 seconds, for example:
    echo "Still running - $(date)"
done
```

find all .gro files and do something with them

```
files=$(ls *gro)
for file in $files
do

#what you do actually want to do 
echo $file
done
```
4) Variable definitions are as simple as e.g. `myvariable=5`. It can be called then by its name adding a $ before. `$myvariable`.

5) Comments: \# is a comment symbol and everything behind it will not be executed.

6) Simple addition `newvalue=$(($oldvalue + $newvalue))`

7) Logic for comparison

| Operator / Command | Meaning / Example |
|---|---|
| `-gt` | Numeric greater than. Example: `[ "$a" -gt "$b" ]` |
| `-lt` | Numeric less than. Example: `[ "$a" -lt "$b" ]` |
| `-eq` | Numeric equal. Example: `[ "$a" -eq "$b" ]` |
| `-ne` | Numeric not equal. Example: `[ "$a" -ne "$b" ]` |
| `-ge` | Numeric greater than or equal. Example: `[ "$a" -ge "$b" ]` |
| `-le` | Numeric less than or equal. Example: `[ "$a" -le "$b" ]` |
| `=` | String equal (POSIX `[` ). Example: `[ "$str1" = "$str2" ]` |
| `!=` | String not equal. Example: `[ "$str1" != "$str2" ]` |
| `<` / `>` | String less/greater (lexicographic) in `[[ ]]` or use `\<'`/`\>'` in `[ ]`. Example: `[[ "$a" < "$b" ]]` |
| `-z` | String is empty. Example: `[ -z "$var" ]` |
| `-n` | String is not empty. Example: `[ -n "$var" ]` |
| `!` | Logical NOT. Example: `if ! [ -f file ]; then ...` |
| `&&` | Logical AND between commands/tests. Example: `[ "$a" -gt 0 ] && echo "pos"` |
| `\|\|` | Logical OR between commands/tests. Example: `[ -f f ] \|\| echo "missing"` |
 | `-e` | File exists (any type). Example: `[ -e path ]` |
| `-f` | Regular file exists. Example: `[ -f file ]` |
| `-d` | Directory exists. Example: `[ -d dir ]` |
| `-r` | File is readable. Example: `[ -r file ]` |
| `-w` | File is writable. Example: `[ -w file ]` |
| `-x` | File is executable. Example: `[ -x script.sh ]` |
| `-s` | File exists and has non-zero size. Example: `[ -s file ]` | 

Notes:
- Use `[ ... ]` (POSIX) or `[[ ... ]]` (bash). `[[ ]]` is more flexible (allows `==`, `<`, `>` without escaping and pattern matching).
- Always quote variables in tests (e.g. `[ "$a" = "$b" ]`) to avoid word-splitting and errors with empty values.
- Numeric comparisons require the `-` operators (`-eq`, `-gt`, …). Using `=`/`!=` compares strings, not numbers.
- For floating-point comparisons use external tools (e.g. `awk`, `bc`) because `test`/`[` only handle integers.
- Combine tests with `&&`, `||`, or use compound tests: `if [ -f a ] && [ -w a ]; then … fi`.
- Use `(( ))` for arithmetic expressions: `if (( a > b )); then … fi` (no `$` needed inside).
- Example: `if (( a > b && b >= 0 )); then echo "ok"; fi`
- Example combining file checks: `if [ -f file ] && [ -s file ]; then echo "exists and not empty"; fi`
8) permissions:
you may somethimes need to change permissions. To do this you can use chmod to change permissions for user, group and the rest. Using a common numbering scheme this means, e.g. `chmod 755` gives you all permissions, whereas everyone else may read and run it but not modify. (common for scripts)

| Permission | Number | chmod shorthand |
|---|---:|---|
| `---` (no permissions) | 0 | `---` |
| `--x` (execute only) | 1 | `--x` |
| `-w-` (write only) | 2 | `-w-` |
| `-wx` (write & execute) | 3 | `-wx` |
| `r--` (read only) | 4 | `r--` |
| `r-x` (read & execute) | 5 | `r-x` |
| `rw-` (read & write) | 6 | `rw-` |
| `rwx` (read, write & execute) | 7 | `rwx` |


## Setting up a Simulation for the combined 1J46.pdb

An MD simulation in GROMACS will always need 3 things
  1) coordinates, which tell the program where is what (usually a `.gro` file)
  2) a "topology" (yes it is defined different in Maths...). This contains the actual forcefield that tells 
  the program how which molecule interacts with all the others and with itself. This can be nicely separated out in a pure `.itp` 
  file for conciseness and to transport it over.
  Furthermore it wants to know (at the beginning) how the atomtypes are defined and how to deal with LR-interactions in the intra-molecular part.
  Also at the end you specify how many of which kind of molecules (in the order of the geometry `.gro`) is included.
  All this makes a `.top`file.

  3) Instructions what you expect the program to actually do with those things. 
  (Energy minimization/ MD run, how to deal with temparature (thermostat parameters), pressure (barostat parameters), constraints,
  how do describe LR-interactions, free energy dH/dlambda calculations,...  )

  Once you have all those, you can send the instructions to the "preprocessor" 
  ``` gmx grompp -f instuctions.mdp -p topology_file.top -c starting_coordinates.gro -o name_how_runfile_should_be_called.tpr```

  This then checks everything for correctness and yields a `.tpr`file.

  This you can then use for a simulation. To simplify the naming conventions and to prevent you from typing differnt stuff for
  energy, trajectory, log, etc. files you can use the `deffnm` option here, which uses the same name for everything and just changes the file endings.
  e.g. this starts a (verbose) run with the name chosen above. (HINT DO NOT WRITE .tpr here, that conflicts with other file endings.)
   ``` gmx mdrun -v -deffnm  name_how_runfile_should_be_called```

  The following are some common GROMACS mdrun output files

  | File ending | Description |
  |---|---|
  | `.gro` | Final finished output geometry: coordinates (and box vectors) in GROMACS `.gro` format, used as the new starting structure. |
  | `.xtc` | Compressed trajectory (coordinates only, reduced precision) for visualization and analysis; small file size. |
  | `.trr` | Full-precision trajectory containing coordinates, velocities, forces; much larger. |
  | `.edr` | Binary energy file with time series of energies and thermodynamic observables; read with `gmx energy -f `. |
  | `.log` | Text log of the mdrun run: settings, progress, performance summary, warnings and errors. |
  | `.cpt` | Checkpoint file for restarting runs: binary snapshot of simulation state (step, random seeds, integrator state); use with `-cpi` to continue. |

  
  
  ### pdb preparation

First get the `1J46.pdb` file from the protein database or our course.

You can quickly have a look both at the visual structure (using `vmd 1J46.pdb`) and the actual text file (e.g. using `more`, `less` or an editor like `vi`, `vim` or `nano`). 

Think of what protonation state do you expect at neutral conditions? (pH=7.2)  What are the amino-acid pKs Values for this? - Does this correspond to your assumptions?

You can use `https://server.poissonboltzmann.org/pdb2pqr` to get an automated generation of charge states of each individual pH for the amino acids. 
Still you have to be careful and you definitely should check its results.

First you want to split the (mixed) pdb file, which contains 2 DNA parts and a protein part.
For this you only need the respective `ATOM` parts in their individual files. 
To speed this up you can use the privided `split_pdb.py` script. This is a very simple python script that iterates over the lines and separates at the separator.
Alternatively you can of course also use a texteditor.

You should be left with 3 `.pdb` files. (Named `part_0.pdb`, `part_1.pdb`, `part_2.pdb` if you used the script to do this) 


### TODO HERE ACTUAL INSTRUCTIONS 

### so far only notes for myself

`make_structure_for_gromacs.sh `so far generates the topology (CHARGES NOT YET SET TO ZERO- DEPENDING ON WHAT WE WANT TO SIMULATE EITHER SIMULATE AT PH, 
chargestate etc)

todo: once finished make this as instructions, so the students do this in part by line and others they can execute


`run_mds.sh`is a script that they will have to do themselves- not for the course but for our usage

## Course graining (day 2) 

## Data evaluation (day 3)
