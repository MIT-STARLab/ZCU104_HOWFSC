rootfolder = system.arguments[0]
solution = system.arguments[1]
logfile = new Stream(rootfolder+'/'+solution+'/'+solution+'.log');
while (!logfile.eof)
{
	l = logfile.readLine()
	if (l.indexOf('Fmax') !== -1) writeln(l)
}
logfile.close();
synth = new Stream(rootfolder+'/'+solution+'/syn/report/csynth.rpt');
while (!synth.eof)
{
	l = synth.readLine()
	if (l.indexOf('Performance & Resource Estimates:') !== -1)
	{
	  for (i=0; i<16; i++) writeln(synth.readLine().trim());
	}
}
synth.close();
for (var i=2; i<system.arguments.length; i++)
{
	kernel = system.arguments[i];
	writeln(kernel)
	cosim = new Stream(rootfolder+'/'+solution+'/sim/report/'+kernel+'_cosim.rpt');
	for (i=0;i<4;i++) cosim.readln();
	while (!cosim.eof) writeln(cosim.readln());
	cosim.close();
}
// type mv_kernel\solution1\sim\report\matrixvector_cosim.rpt 