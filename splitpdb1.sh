### Author: Michael Dolan
### Date:  2003 ???
### Usage: split_pdb.sh <multi-pdb file> [<file prefix>]     
###
###        <multi-pdb file> Name of the multi-pdb file to split
###        <file prefix>     File prefix for split files (default file)
###
### Description:
###
### Splits a PDB file into separate files
###

#!/bin/sh

usage='usage: split_pdb.sh <multi-pdb file> [<file prefix>]'

if [ $# -lt 1 -o $# -gt 2 ];then
   echo $usage
   exit 1
fi

if [ ! -r $1 ];then
   echo "Multi-pdb file $1 not found or is empty"
   exit 2
fi

if [ $# -eq 2 ];then
   file=$2
else
   file="file"
fi

#echo "using prefix: ${file}"
# nawk can be used here as well

awk '
BEGIN 	{ 
   n =  0
   fnam = sprintf( "%01d.pdb", n=n )
   at_top  = "yes"
   comment = "no"
}
/^#/  	{
   if ( comment == "yes" ) { 
      at_top = "no"
   } else {
      at_top  = "yes"
   }
   comment = "yes"
}
/END/ {
   if ( comment == "no" ) { 
      at_top = "yes"
   } else {
      at_top = "no"
   }
   comment = "no"
}
{ 
   if ( at_top == "yes" ) { 
      close( fnam )
      fnam = sprintf( "%01d.pdb", n=n+1 )
      at_top = "no"
   }
   print $0 > fnam
}' $1

#perl -pi -e 's/END/#/g' *.pdb
