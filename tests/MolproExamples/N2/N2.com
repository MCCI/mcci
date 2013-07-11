 MEMORY, 250, M
 R1=2.0
 geometry={N;N,N,R1}
 basis=6-31G
 hf;
{FCI;DUMP}
