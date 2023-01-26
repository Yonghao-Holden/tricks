import sys;



ifpath=sys.argv[1];
ofpath=sys.argv[2];
ofile=open(ofpath, "w+")

nIndex=0;
for line in open( ifpath ):
	nIndex=nIndex+1;
	if nIndex%4 != 1:
		continue;

	strArr=line.split(" ");
	strName=strArr[0];
	strExtra=strArr[1];
	strExtras=strExtra.split(";");

	strUMI=strExtras[1].lstrip("FP:");

	ofile.write( strName+"\n");
	ofile.write( strUMI+"\n");
	ofile.write( "+\n");
	ofile.write( "!"*len(strUMI)+"\n");

	

#	@M02084:463:000000000-CCKBV:1:1101:17735:1465 1:N:0:GTTATCCCT;FP:CCACGGGTG;RQ:32.59

ofile.close();



 
