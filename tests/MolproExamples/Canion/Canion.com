***,File
 memory,200,m
 geometry={C
 }
 set,charge=-1
 set,spin=3
 basis=vdz
 hf;
 {FCI;CORE,1;DUMP}
