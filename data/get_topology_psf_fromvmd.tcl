#
proc newRep { sel type color rep imol} {
  
    mol selection $sel
    mol representation $type
    mol addrep $imol
    mol showrep $imol $rep on
    mol modcolor $rep $imol $color

}


set PDBFILE "./pe_1mon_aa.pdb"
mol new $PDBFILE type pdb
set imol1 [molinfo top]
mol delrep 0 $imol1
# ALL ATOMS
set rep 0
newRep "all" "CPK " "Name" $rep $imol1

topo getanglelist

animate write psf pe_1mon_aa.psf




