proc newRep { sel type color rep imol} {
  
    mol selection $sel
    mol representation $type
    mol addrep $imol
    mol showrep $imol $rep on
    mol modcolor $rep $imol $color

}

display projection orthographic
axes location off
color Display Background black
display depthcue off

# Open file containing the ligands and read it line by line
set f [open "ligands_list.txt" r]
set lines [split [read $f] "\n"]
close $f
# Read line by line
foreach line $lines {
  set pdb [lindex [split $line " "] 0]
  set reslig [lindex [split $line " "] 1]
  
  mol new $pdb type {webpdb} first 0 last -1 step 1 waitfor 1
  set imol1 [molinfo top]
  mol delrep 0 $imol1
  set rep 0
  if {$reslig == ""} {
	  newRep "not protein and not water and conformationA" "CPK" "element" $rep $imol1
	  set sel [atomselect top "not protein and not water and conformationA"]
  } else {
	  newRep "resid $reslig" "CPK" "ColorID 0" $rep $imol1
	  set sel [atomselect top "resid $reslig"]
  }
  $sel writepdb ${pdb}_ligand.pdb
  mol delete top

}

