#! /usr/bin/perl 
# DMR: 16/12/08.
# NOTE: Place propka2.0, dgl.py and pIl.py in your working directory (pwd) or in your bin directory (~/bin). OPEN BABEL is also required.
# NOTE: Alphabetic alternate atom labels are relabelled as integers in this script.

Usage() unless (@ARGV); 

# Look for propka2.0, dg.py and pI.py, and babel
 $lok=`which propka2.0`;chomp($lok);
 if((stat($lok)) || (stat("propka2.0"))) {;} else {die ("Couldn't find propka2.0")};
 $lok=`which dgl.py`;chomp($lok);
 if((stat($lok)) || (stat("dgl.py"))) {;} else {die ("Couldn't find dgl.py")};
 $lok=`which pIl.py`;chomp($lok);
 if((stat($lok)) || (stat("pIl.py"))) {;} else {die ("Couldn't find pIl.py")};
 $lok=`which babel`;chomp($lok);
 if((stat($lok)) || (stat("babel"))) {;} else {die ("Couldn't find babel")};

# Read input
for($inp=0;$inp<=$#ARGV;$inp++){
 if($ARGV[$inp] eq '-i'){$input=$ARGV[$inp+1];}
 if($ARGV[$inp] eq '-c'){$chain=$ARGV[$inp+1];}
 if($ARGV[$inp] eq '-m'){$domodi=$ARGV[$inp+1];}
 if($ARGV[$inp] eq '-s'){$salt=1;}
 if($ARGV[$inp] eq '-apo'){$apo=1;}
                               }
 if($ARGV[0] ne '-i'){print "Please input -i NAME.pdb\n"; exit;}

 if($chain ne ''){
 @chains=split(/,/,$chain);
 if($#chains == 0){print "1 chain ID specified: @chains \n"};
 if($#chains > 0){print "More than 1 chain ID specified: @chains \n"};
                   } $chain=~s/\,//g;

# Read PDB file
open (file,"$input") or die ("Couldn't open file $input, stopped");
$ic=0;
while(<file>) {
 $ic++;
 chomp;
 $atom[$ic]=$_;
 if(substr($atom[$ic],0,4)=~/^ATOM/ or substr($atom[$ic],0,6)=~/^HETATM/ and substr($atom[$ic],16,1)=~/\w/){substr($atom[$ic],16,1)=~tr/A-Z/1-26/;} #Label alternate atom A as atom 1,... 
              }
 if($atom[$ic]=~/^LG/){rprop($input);} #Call sub rprop if Propka2.0 input file
# Process PDB file

                  $code=substr($input,0,-4); $len=$#atom; $ic=0; $imod=0; $id=0; $ir=0;
	      for($j=1;$j<=$len;$j++){
	      if($atom[$j]=~/^SEQRES/ and substr($atom[$j],8,2)=~/\s1/){$ic++;$chn[$ic]=substr($atom[$j],11,1);} #CHAIN ID column 12
	      if($atom[$j]=~/^MODEL/ and $imod==0){modstr(@atom,$j,$imod,@mod,@mods,@mode);} #Call subroutine modstr if model structure(s) in PDB file
		                     }
	      if($imod>=1){$mlen=$#mod;print "$imod Models/Structures in PDB file $input. \n";}
	      if($imod==0){$mlen=0;}


	      if($imod==0 or $imod>=1){for($im=0;$im<=$mlen;$im++){$domod=$mod[$im]; #Do model loop
              if($domodi=~/\w/){$domod=$domodi;$im=$mlen+1;} #Do input model $domodi

	      if($imod==0){$len=$#atom;$iddup=0;$str=1;$end=$len;}
	      if($imod>=1){$len=$#atom;$iddup=0;$str=$mods[$domod-1];$end=$mode[$domod-1];
	         print "Doing Model $domod; start $str, end $end.\n";}


	      for($j=$str;$j<=$end;$j++){

	      if($atom[$j]=~/^ATOM/ and $ic==0){nohead(@atom,$j,$ic,@chn)}; #Call subroutine nohead if NO SEQRES info in PDB file

	      if($atom[$j]=~/^ATOM/ and substr($atom[$j],16,1)!~/\s/){$id++;$dup[$id]=substr($atom[$j],16,1);$dup2[$id]=substr($atom[$j],22,4);if($dup2[$id]!=$dup2[$id-1]){$dupres++;push(@dup3,$dup2[$id])}} #Check for duplicate atoms in column 17

	      if($atom[$j]=~/^ATOM/ and substr($atom[$j],26,1)!~/\s/ and substr($atom[$j],22,4)==substr($atom[$j-1],22,4) and substr($atom[$j],22,5)!~substr($atom[$j-1],22,5)){$ir++;$reir[$ir]=substr($atom[$j-1],22,4);} #Check for residues with same resSeq no.

	                                }

# Sort duplicate atom labels
	      @dups = sort $a <=> $b, @dup; 
	      for($i=1;$i<=$#dups;$i++){if($dups[$i+1]!~$dups[$i]){$iddup++;push(@conf,$dups[$i]); }} @confp=@conf;shift(@confp);

# Status printing
	      if($id>=1){print "WARNING: ".@confp." alternate conformer sets found; "."@confp"." for ".$id." atoms in ".$dupres." residues,"."@dup3"."\nNOTE: Alphabetic alternate atom labels are relabelled as integers.\nNOTE: Only the 1st conformer set will be used.\n";}
	      if($ir>=1){print "WARNING: Residues with same resSeq number,"."@reir".". Please check PDB file.\n";}
	      if($ic==1){print "One ATOM chain ".$ic.", "."@chn".". Doing chain ".$chain."\n";}
	      if($ic>1){print "Number of ATOM chains ".$ic.", "."@chn"; if($chain=~/\w/){print  ". Doing chain ".$chain."\n";} if($chain!~/\w/){print". Doing all\n";}}


              if($chain eq ''){$ic=1};
                 $clim=$#conf;


# Generate temporary/new PDB files

	      if($dupres==0){ #No duplicate atoms

	      if($ic==1){ #1 chain
              $output="$code\_all.pdb";open (file2,">$output");
              $ja=0;
              for($j=$str;$j<=$end;$j++){
              if($atom[$j]=~/^ATOM/                        ){ 
              if(substr($atom[$j],13,1)=~/H/){next;} #skip H column 14
              if(substr($atom[$j],12,1)=~/H/){next;} #skip H column 13
               $ja++;$atmp[$ja]=$atom[$j];
        	      print file2 $atmp[$ja]."\n";
              } } 
if($apo!=1){
              for($j=$str;$j<=$end;$j++){
              if(                      $atom[$j]=~/^HETATM/){ 
              if(substr($atom[$j],13,1)=~/H/){next;} #skip H column 14
              if(substr($atom[$j],12,1)=~/H/){next;} #skip H column 13
              if(substr($atom[$j],17,3)=~/HOH/){next;} #skip HOH 
               $ja++;$atmp[$ja]=$atom[$j];
        	      print file2 $atmp[$ja]."\n";
	      } }  }
			}#1 chain


	      if($#chains>=0){ #do input chains
              $output="$code\_$chain.pdb";open (file2,">$output");

 	      while(<@chains>){$chan=$_;print"$chan\n";

              $ja=0;
              for($j=$str;$j<=$end;$j++){
              if($atom[$j]=~/^ATOM/                        ){ 
              if(substr($atom[$j],21,1)!~/$chan/ and $atom[$j]=~/^ATOM/){next;} #Chain ID col 22
              if(substr($atom[$j],13,1)=~/H/){next;} #skip H column 14
              if(substr($atom[$j],12,1)=~/H/){next;} #skip H column 13
              $ja++;$atmp[$ja]=$atom[$j];
        	      print file2 $atmp[$ja]."\n";
              } } 
			    }#chains
if($apo!=1){
 	      while(<@chains>){$chan=$_;#print"$chan\n"; #added

              for($j=$str;$j<=$end;$j++){
              if(                      $atom[$j]=~/^HETATM/){ 
              if(substr($atom[$j],21,1)!~/$chan/){next;} #skip 
              if(substr($atom[$j],13,1)=~/H/){next;} #skip H column 14
              if(substr($atom[$j],12,1)=~/H/){next;} #skip H column 13
              if(substr($atom[$j],17,3)=~/HOH/){next;} #skip HOH 
               $ja++;$atmp[$ja]=$atom[$j];
        	      print file2 $atmp[$ja]."\n";
	      } }
			    }  }#chains
			    }#do input chains

			    }#No duplicates
			    


	      if($dupres!=0){ #Duplicate atoms

	      if($ic==1){ #1 chain
              for($icf=1;$icf<=$clim;$icf++){ 
              $confp=$conf[$icf];print "Doing alternate atom set $confp\n"; 
              $output="$code\_conf$confp.pdb";open (file2,">$output");

              $ja=0; 
              for($j=$str;$j<=$end;$j++){
	      $iop=0;
              if($atom[$j]=~/^ATOM/                        ){ 
              if(substr($atom[$j],13,1)=~/H/){next;} #skip H column 14
              if(substr($atom[$j],12,1)=~/H/){next;} #skip H column 13
              if($iop!=1 and substr($atom[$j],16,1)!~/\s/ and substr($atom[$j],16,1)ne$conf[$icf]){next;}#skip if altLoc ne conf
              $ja++;$atmp[$ja]=$atom[$j];
        	      print file2 $atmp[$ja]."\n";
              } } 
if($apo!=1){
              for($j=$str;$j<=$end;$j++){
	      $iop=0;
              if(                      $atom[$j]=~/^HETATM/){ 
              if(substr($atom[$j],13,1)=~/H/){next;} #skip H column 14
              if(substr($atom[$j],12,1)=~/H/){next;} #skip H column 13
              if(substr($atom[$j],17,3)=~/HOH/){next;} #skip HOH 
              if($iop!=1 and substr($atom[$j],16,1)!~/\s/ and substr($atom[$j],16,1)ne$conf[$icf]){next;}#skip if altLoc ne conf
               $ja++;$atmp[$ja]=$atom[$j];
        	      print file2 $atmp[$ja]."\n";
	      } }  }

				             }#icf
			}# 1 chain


	      if($#chains>=0){ #do input chains
              for($icf=1;$icf<=$clim;$icf++){ 
              $confp=$conf[$icf];print "Doing alternate atom set $confp\n";
              $output="$code\_$chain\_conf$confp.pdb";open (file2,">$output");

 	      while(<@chains>){$chan=$_;print"$chan\n";

              $ja=0;
              for($j=$str;$j<=$end;$j++){
	      $iop=0;
              if($atom[$j]=~/^ATOM/                        ){ 
              if(substr($atom[$j],21,1)!~/$chan/ and $atom[$j]=~/^ATOM/){next;} #Chain ID col 22
              if(substr($atom[$j],13,1)=~/H/){next;} #skip H column 14
              if(substr($atom[$j],12,1)=~/H/){next;} #skip H column 13
              if($iop!=1 and substr($atom[$j],16,1)!~/\s/ and substr($atom[$j],16,1)ne$conf[$icf]){next;}#skip if altLoc ne conf
              $ja++;$atmp[$ja]=$atom[$j];
        	      print file2 $atmp[$ja]."\n";
              } } 
           			     }#chains
if($apo!=1){
 	      while(<@chains>){$chan=$_;#print"$chan\n"; #added

              for($j=$str;$j<=$end;$j++){
	      $iop=0;
              if(                      $atom[$j]=~/^HETATM/){ 
              if(substr($atom[$j],21,1)!~/$chan/){next;} #skip 
              if(substr($atom[$j],13,1)=~/H/){next;} #skip H column 14
              if(substr($atom[$j],12,1)=~/H/){next;} #skip H column 13
              if(substr($atom[$j],17,3)=~/HOH/){next;} #skip HOH 
              if($iop!=1 and substr($atom[$j],16,1)!~/\s/ and substr($atom[$j],16,1)ne$conf[$icf]){next;}#skip if altLoc ne conf
               $ja++;$atmp[$ja]=$atom[$j];
        	      print file2 $atmp[$ja]."\n";
	      } }

           			     }  }#chains
	          			     }#icf
			}#do input chains

			    }#duplicate atoms

close(file1,file2);

# Run PROPKA and Python scripts

	      if($dupres==0){ #No duplicate atoms
# Build new_PDB.pdb file using $output from above (as DCB's pre_propka.sh)
if($apo!=1){ prepropka($output); $output="new_$output"; }
# Add model pKa ligand values
if($apo!=1){ spkamod($output,$salt); }
#Run propka (in ~/bin)
 $JOB=$output; print $JOB."\n";
#system("rm $JOB.pka");
 system("cp $JOB PROPKATMP2");
 system("propka2.0");
#system("rm PROPKATMP1 PROPKATMP2 $output");
 system("rm            PROPKATMP2");
 if($imod==0){
 if($ic==1){ system("cp PROPKATMP3 $code.pdb.pka");$pkaout="$code.pdb.pka";}
 if($ic>=2){ system("cp PROPKATMP3 $code\_$chain.pdb.pka");$pkaout="$code\_$chain.pdb.pka";}
 }
 if($imod>=1){
 if(($im==0 or $im==$mlen+1) and $ic==1){ system("echo 'Model $domod' > $code.pdb.pka && cat PROPKATMP3 >> $code.pdb.pka                 ");$pkaout="$code.pdb.pka";}
 if($im>=1 and $im<=$mlen and $ic==1){ system("echo '\nModel $domod' >> $code.pdb.pka && cat PROPKATMP3 >> $code.pdb.pka                 ");$pkaout="$code.pdb.pka";}
 if(($im==0 or $im==$mlen+1) and $ic>=2){ system("echo 'Model $domod' > $code\_$chain.pdb.pka && cat PROPKATMP3 >> $code\_$chain.pdb.pka                 ");$pkaout="$code\_$chain.pdb.pka";}
 if($im>=1 and $im<=$mlen and $ic>=2){ system("echo '\nModel $domod' >> $code\_$chain.pdb.pka && cat PROPKATMP3 >> $code\_$chain.pdb.pka                 ");$pkaout="$code\_$chain.pdb.pka";}
 }

#Run dg.py and pI.py (in ~/bin)
 open(file3,"PROPKATMP3");open(file4,">tmp.pka");
 while(<file3>) {
 chomp;
 if($_=~/^\sWarning/){system("rm PROPKATMP3 tmp.pka"); die ("Calculation of free energy and pI stopped: Missing ionizable residue atoms (see PROPKA 2.0 output).\n");}
 if($_=~/^\s\s\s\w/ and $_!~/^\s\s\sRESIDUE/ and substr($_,13,5)!~/99.99/){print file4 $_."\n";}
                }

 open(file3,">>$pkaout");open(file4,"tmp.pka");
 $pyout=`dgl.py tmp.pka`; print $pyout; print file3 $pyout; $pyout=`pIl.py tmp.pka`; print $pyout; print file3 $pyout;
 
 close(file3,file4); system("rm PROPKATMP3 tmp.pka");
			    }


	      if($dupres!=0){ #Duplicate atoms
              for($icf=1;$icf<=$clim;$icf++){ 
		  $confp=$conf[$icf];
              if($ic==1){$output="$code\_conf$confp.pdb";open (file2,"$output");}
              if($ic>=2){$output="$code\_$chain\_conf$confp.pdb";open (file2,"$output");}
# Build new_PDB.pdb file using $output from above (as DCB's pre_propka.sh)
if($apo!=1){ prepropka($output); $output="new_$output"; }
# Add model pKa ligand values
if($apo!=1){ spkamod($output,$salt); }
#Run propka (in ~/bin)
 $JOB=$output; print $JOB."\n";
#system("rm $JOB.pka");
 system("cp $JOB PROPKATMP2");
 system("propka2.0");
#system("rm PROPKATMP1 PROPKATMP2 $output");
 system("rm            PROPKATMP2");
 if($imod==0){
 if($ic==1){ system("cp PROPKATMP3 $code\_conf$confp.pdb.pka");$pkaout="$code\_conf$confp.pdb.pka";}
 if($ic>=2){ system("cp PROPKATMP3 $code\_$chain\_conf$confp.pdb.pka");$pkaout="$code\_$chain\_conf$confp.pdb.pka";}
 }
 if($imod>=1){
 if(($im==0 or $im==$mlen+1) and $ic==1){ system("echo 'Model $domod' > $code\_conf$confp.pdb.pka && cat PROPKATMP3 >> $code\_conf$confp.pdb.pka                 ");$pkaout="$code\_conf$confp.pdb.pka";}
 if($im>=1 and $im<=$mlen and $ic==1){ system("echo '\nModel $domod' >> $code\_conf$confp.pdb.pka && cat PROPKATMP3 >> $code\_conf$confp.pdb.pka                 ");$pkaout="$code\_conf$confp.pdb.pka";}
 if(($im==0 or $im==$mlen+1) and $ic>=2){ system("echo 'Model $domod' > $code\_$chain\_conf$confp.pdb.pka && cat PROPKATMP3 >> $code\_$chain\_conf$confp.pdb.pka                 ");$pkaout="$code\_$chain\_conf$confp.pdb.pka";}
 if($im>=1 and $im<=$mlen and $ic>=2){ system("echo '\nModel $domod' >> $code\_$chain\_conf$confp.pdb.pka && cat PROPKATMP3 >> $code\_$chain\_conf$confp.pdb.pka                 ");$pkaout="$code\_$chain\_conf$confp.pdb.pka";}
 }

#Run dg.py and pI.py (in ~/bin)
 open(file3,"PROPKATMP3");open(file4,">tmp.pka");
 while(<file3>) {
 chomp;
 if($_=~/^\sWarning/){system("rm PROPKATMP3 tmp.pka"); die ("Calculation of free energy and pI stopped: Missing ionizable residue atoms (see PROPKA 2.0 output).\n");}
 if($_=~/^\s\s\s\w/ and $_!~/^\s\s\sRESIDUE/ and substr($_,13,5)!~/99.99/){print file4 $_."\n";}
                }
		                         
 open(file3,">>$pkaout");open(file4,"tmp.pka");
 $pyout=`dgl.py tmp.pka`; print $pyout; print file3 $pyout; $pyout=`pIl.py tmp.pka`; print $pyout; print file3 $pyout;

 close(file3,file4); system("rm PROPKATMP3 tmp.pka");
                                             }
			    } 

			    }} #Do model loop

 exit;


              sub Usage() {
  print "propka2.0.pl -i <PDB input file> -c <protein and ligand chain IDs separated by ',' (OR skip/leave blank for all chains)> -s (include salts as charged atoms) -apo (exclude all HETATMs) -m <NMR Model (if applicable)> \n";
  print "Note: propka2.0.pl assumes protein is ATOM and ligand is HETATM record type; H atoms, HOH are ignored.\n";
  print "Examples:\n";
  print "propka2.0.pl -i 1K1I.pdb (chain A w/ ligand FD1)\n";
  print "propka2.0.pl -i 4DFR.pdb -c B (chain B w/ ligand MTX)\n";
  print "propka2.0.pl -i 4ER2.pdb -c E,I (chain E w/ ligand chain I, pepstatin. Use edited PDB file.)\n"; exit;
    }


	      sub nohead(@atom,$j,$ic,$chn){
	      $len=$#atom;$ic=0;
	      for($i=$j;$i<=$len;$i++){
      if($atom[$i]=~/^ATOM/){ push @chnh, substr($atom[$i],21,1);}
    }
              for($i=0;$i<=$#chnh;$i++){
      if($chnh[$i+1]!~$chnh[$i]){$ic++;push(@chn,$chnh[$i]);}
    } 
    }


	      sub modstr(@atom,$j,$imod,@mod,@mods,@mode){
	      $len=$#atom;$imod=0;
	      for($i=$j;$i<=$len;$i++){
      if($atom[$i]=~/^MODEL/){ $imod++; push @mod, substr($atom[$i],10,4);}
      if($atom[$i]=~/^MODEL/){ push @mods, $i;}
      if($atom[$i]=~/^ENDMDL/){ push @mode, $i;}
    }
    }


    sub prepropka($output){
#From DCB's pre_propka.sh. Re-written in script, incl. conv_type.f (below).
 my @tmp1,@tmp2;
 $JOB=$output;
#
#
# CREATION OF PDB FILE FOR THE LIGAND ONLY
#
#
 open(filei,"$output"); open(fileo,">new_ligand.pdb"); open(filet,">TMP1");
 while(<filei>){
  if($_ =~ /^HETATM/){
	 	 print fileo $_;
	 if(substr($_,12,1)=~/\w/){substr($_,13,2)=substr($_,12,2);substr($_,12,1)=' ';}
                 print filet $_;
 } } close(filei,fileo,filet);
#
#
# TWO-STEPS CONVERSION TO THE SYBYL2 FORMAT
#
#
 system("babel new_ligand.pdb new_ligand.xyz"); #ok using open babel 2.1.1
 system("babel new_ligand.xyz new_ligand.mol2"); 
#   
 open(fileli,"new_ligand.mol2"); open(filelo,">TMP2");
 while(<fileli>){
 chomp, $_;
 @tmp = split, $_;
 if($tmp[0]=~/[1-9]/ and $tmp[1]=~/[A-Z]/ and $tmp[1]!~/^H/ and $tmp[1]!~/^X/) {$tmp[5]=~s/\.//g; print filelo $tmp[5],"\n";}
 }
 close(fileli,filelo);
#
#
 open(filepi,"$JOB"); open(filepo,">TMP3");
 while(<filepi>){
 chomp, $_;
 if($_=~/^ATOM/) {print filepo $_,"\n";}
 }
 close(filepi,filepo);
#
#DMR: Below replaces call to DCB's conv_type.f
 open(filei,"TMP3"); open(fileo,">protein_new");
 while(<filei>){
 chomp, $_;
 substr($_,60,6)=' 00.00'; 
 printf fileo "%-81s\n",$_;
 }
 close(filei,fileo);
 open(filei1,"TMP1"); open(filei2,"TMP2"); open(fileo,">PDBNEW");
 @tmp1=();@tmp2=();
 while(<filei1>){
 chomp, $_;
 substr($_,60,6)=' 00.00'; 
 push @tmp1, $_;
 }
 while(<filei2>){
 chomp, $_;
 push @tmp2, $_;
 }
 for($i=0;$i<=$#tmp1;$i++){
 printf fileo "%2s%4s%7s%-3s %-64s\n",'LG',substr($tmp1[$i],13,3),substr($tmp1[$i],6,7),substr($tmp2[$i],0,3),substr($tmp1[$i],17);
 }
 close(filei1,filei2,fileo);
# 
# CREATION OF THE PDB FILE FOR THE PROTEIN
#
 system("cat protein_new PDBNEW > new_$JOB");
 system("rm TMP*");
 system("rm new_ligand.*");
 system("rm PDBNEW");
 system("rm protein_new");
 }


 sub spkamod($output,$salt){
#Substitute ligand pKas with pKa model values
 $filei=$output;
 open(filei,"$filei"); open(fileo,">spkamod.tmp");
 while(<filei>){
       if(substr($_,0,4)=~/ATOM/){print fileo $_;}
       if(substr($_,0,2)=~/LG/){
	       if(substr($_,13,3)=~/Npl/){substr($_,60,6)=' 11.50';} 
	       if(substr($_,13,3)=~/N3 /){substr($_,60,6)=' 10.00';} 
	       if(substr($_,13,3)=~/N4 /){substr($_,60,6)=' 10.00';} #830C.pdb
	       if(substr($_,13,3)=~/Nar/){substr($_,60,6)=' 05.00';} 
	       if(substr($_,13,3)=~/Oco/){substr($_,60,6)=' 04.50';} 
	         if($salt){
      #if(substr($_,13,3)=~/Li / and substr($_,17,3)=~/ LI/){substr($_,13,3)='1P ';} #Li 1+
      #if(substr($_,13,3)=~/Na / and substr($_,17,3)=~/ NA/){substr($_,13,3)='1P ';} #Na 1+
       if(substr($_,13,3)=~/Mg / and substr($_,17,3)=~/ MG/){substr($_,13,3)='2P ';} #Mg 2+
       if(substr($_,13,3)=~/K  / and substr($_,17,3)=~/  K/){substr($_,13,3)='1P ';} #K  1+
       if(substr($_,13,3)=~/Ca / and substr($_,17,3)=~/ CA/){substr($_,13,3)='2P ';} #Ca 2+
       if(substr($_,13,3)=~/Fe / and substr($_,17,3)=~/ FE/){substr($_,13,3)='2P ';} #Fe 2+
       if(substr($_,13,3)=~/Ni / and substr($_,17,3)=~/ NI/){substr($_,13,3)='2P ';} #Ni 2+
       if(substr($_,13,3)=~/Cu / and substr($_,17,3)=~/ CU/){substr($_,13,3)='2P ';} #Cu 2+
       if(substr($_,13,3)=~/Zn / and substr($_,17,3)=~/ ZN/){substr($_,13,3)='2P ';} #Zn 2+
       if(substr($_,13,3)=~/Mn / and substr($_,17,3)=~/ MN/){substr($_,13,3)='2P ';} #Mn 2+
      #if(substr($_,13,3)=~/Cl / and substr($_,17,3)=~/ CL/){substr($_,13,3)='1N ';} #Cl 1-
      #if(substr($_,13,3)=~/S3 / and substr($_,17,3)=~/SO4/){substr($_,13,3)='2N ';} #SO4 2-
       }
	       print fileo $_;}
	      }
	      system("mv spkamod.tmp $filei");
      }


 sub rprop($input){
 system("cp $input PROPKATMP2");
 system("propka2.0");
 open(file3,"PROPKATMP3");open(file4,">tmp.pka");
 while(<file3>) {
 chomp;
 if($_=~/^\sWarning/){system("mv PROPKATMP3 $input.pka && rm PROPKATMP2 tmp.pka"); die ("Calculation of free energy and pI stopped: Missing ionizable residue atoms (see PROPKA 2.0 output).\n");}
 if($_=~/^\s\s\s\w/ and $_!~/^\s\s\sRESIDUE/ and substr($_,13,5)!~/99.99/){print file4 $_."\n";}
                }
 open(file3,">>PROPKATMP3");open(file4,"tmp.pka");
 $pyout=`dgl.py tmp.pka`; print $pyout; print file3 $pyout; $pyout=`pIl.py tmp.pka`; print $pyout; print file3 $pyout;
 close(file3,file4);
 system("mv PROPKATMP3 $input.pka && rm PROPKATMP2 tmp.pka");
 exit;
 }
