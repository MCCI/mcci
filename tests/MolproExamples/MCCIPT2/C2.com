 MEMORY, 250, M
 R1=1.5
 geometry={C;C,C,R1}
 basis=6-31G
 hf;
{FCI;CORE;DUMP}
