PROPKA 2.0
==========

PROPKA (http://propka.ki.ku.dk/) 16/12/08

Files included: propka2.0_2008-11-12.f, propka2.0.pl, dgl.py and pIl.py.

propka2.0_2008-11-12.f: PROPKA 2.0 fortran source code.


# propka2.0.pl:
Executes PROPKA 2.0 and runs the python scripts dgl.py and pIl.py. Requires that the PROPKA 2.0 executable file is named 'propka2.0' and located in the pwd or in a bin directory along with dgl.py and pIl.py. In addition [Open Babel](http://www.openbabel.org) is required.

Usage:  propka2.0.pl -i xxxx.pdb [-c 'Chain ID(s)' -s -apo -m 'NMR Model Number']
Example:  propka2.0.pl -i 1K1I.pdb (chain A w/ ligand FD1) 


# dgl.py:
Calculates pH dependent free energy using PROPKA output.

Usage:  dgl.py xxxx.pdb.pka


# pIl.py:
Calculates the pI (isoelectric point) of the protein using PROPKA output.

Usage:  pIl.py xxxx.pdb.pka


# Note:

propka2.0.pl can execute PROPKA 2.0 using a previously generated input file with edited ligand pKmodel values (new_xxxx.pdb).

Usage:  propka2.0.pl -i new_xxxx.pdb
Example:  propka2.0.pl -i new_1K1I.pdb (chain A w/ ligand FD1) 

