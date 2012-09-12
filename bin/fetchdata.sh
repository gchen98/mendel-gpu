#!/bin/bash

if [ $# -lt 4 ] ; then
  echo "Usage <dbname> <chr> <regionsize> <subjectchunksize>"
  exit 1
fi

dbname=$1
chr=$2
regionsize=$3
chunksize=$4

basepath=$IMPUTATION_DATA'/'$dbname
chunks=10

cd $IMPUTATION_DATA
echo Fetching data on $dbname for Chr $chr
mysql -h $DBHOST -u $DBUSER --password=$DBPW $dbname << END |chunk $chr $basepath $chunksize
select count(*) from person;
select chrom,position,rs,vector from $MASTERDB.$MASTERTABLE as a,genotype as b where chrom='$chr' and position>0 and a.rsnumber=b.rs order by position;
END
for((chunk=0;chunk<$chunks;++chunk))
do
  basefile=$basepath'/chr.'$chr'.chunk.'$chunk
  if [ -f $basefile ] ; then
    split_file.pl $basefile'.part' $regionsize 0 < $basefile 
    gzip  -f $basefile'.part'*
    rm $basefile
  fi
done