file=$1

params=amber99

tinkerdir=..
bindir=$tinkerdir/bin
paramdir=$tinkerdir/params
paramfile=$paramdir/$params
pdbhstrip=pdbhstrip
tinkerpatch=./tinkerpatch

basefile=`basename $file .pdb`
basefile=`basename $basefile .ent`
\rm -f $basefile.xyz* $basefile.pdb_* $basefile.seq* foo.*
seqfile=$basefile.seq
xyzfile=$basefile.xyz
xyz2file=${xyzfile}_2
resultfile=${basefile}_result.pdb


pdbxyz=$bindir/pdbxyz
minimize=$bindir/minimize
xyzpdb=$bindir/xyzpdb

# Convert to xyz
### We need to check for alternate occupancies first ###
$pdbxyz $file ALL ALL $paramfile

# Cartesian minimization
$minimize $xyzfile $paramfile 2

# Convert results to pdb
cp $xyz2file foo.xyz
cp $seqfile  foo.seq
$xyzpdb foo.xyz $paramfile
$tinkerpatch $file foo.pdb | $pdbhstrip > $resultfile
rm -f foo.xyz foo.seq foo.pdb

# Remove intermediate files
rm $xyzfile $xyz2file $seqfile

