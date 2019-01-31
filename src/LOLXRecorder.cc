#include "LOLXRecorder.hh"

	LOLXRecorder::LOLXRecorder() : LXeRecorderBase() {
		fout=NULL;
        tree=NULL;
        event=NULL;
        LOLXOutputFileName = "Outputfile.root";
        counter = 0;
	}
