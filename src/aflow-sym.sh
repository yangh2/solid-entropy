#!/bin/bash

usage() { echo "Usage: $0  [ -x <xyz file> ] [ -p <POSCAR file> ]" 1>&2; exit 1;}

pos=1;
fpos="POSCAR";

while getopts ":x:p:" f;
do
    case "$f" in
	x)
	    fxyz=$OPTARG
	    pos=0;
	    ;;
	p)
	    fpos=$OPTARG
	    pos=1;
	    ;;
	*)
	    usage
	    ;;
    esac
done

cwd=`pwd`;

if [ $pos -eq 0 ];
then
    dirname=`mktemp -d`;
    cp $fxyz $dirname;
    cd $dirname;
    xyz2pos $fxyz;
    aflow --aflowSYM < POSCAR 1>tmp;
    xz -d aflow.agroup.out.xz;
    cp aflow.agroup.out $cwd;
    cd $cwd;
    rm $dirname -r;
fi

if [ $pos -eq 1 ];
then
    dirname=`mktemp -d`;
    cp $fpos $dirname;
    cd $dirname;
    aflow --aflowSYM < $fpos 1>tmp;
    xz -d aflow.agroup.out.xz;
    cp aflow.agroup.out $cwd;
    cd $cwd;
    rm $dirname -r;
fi

outfile="symmetry_operation.dat";

if [ -e $outfile ];
then
    rm $outfile;
fi

touch $outfile;
awk 'NR==4{print $1}' aflow.agroup.out >> $outfile;
grep Uf aflow.agroup.out -B 2 | tail -n +3 |sed 's/Uf//g' | awk '$1!="--"{print $0}' >> $outfile;


