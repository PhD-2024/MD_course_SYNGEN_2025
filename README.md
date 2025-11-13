# MD_course

# Linux introduction 

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
| `echo` | prints to stout | 



1) (Bash) Scripts are essentially nothing essentially nothing more storing commands in a text file so you can reuse them (also great to have reproducability - no arbitray behavior due to typos-)
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
While condition fullfilled (file exists)

TODO add 

find all .gro files and do something with them

```
files=$(ls *gro)
for file in $files
do

\#what you do
done
```

## Setting up a Simulation for the combined 1J46.pdb

## pdb preparation

First get the `1J46.pdb` file from the protein database or our course.

You can quickly have a look both at the visual structure (using `vmd 1J46.pdb`) and the actual text file (e.g. using `more`, `less` or an editor like `vi`, `vim` or `nano`)