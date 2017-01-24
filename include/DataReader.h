#ifndef DATAREADER_H
#define DATAREADER_H
/*
 *        Class: DataReader
 *  Description: handling of SAMBA binarary file in streamer mode
 *       Author: Federico Ferri, CEA/Saclay
 */
#include "Types.h"

#include "TFile.h"
#include "TObjString.h"
#include "TTree.h"

class DataReader
{
        public:
                DataReader(const char * output_file_name);

                ~DataReader();

                void readStreamerModeFile(const char * input_file_name);

        private:
                TTree * _t;
                TBranch * _br;
                TFile * _fout;
                FILE *  _fd;
                size_t _nsamples = 1000;
                size_t _ndetids = 0;
                int32_t _detid;
                Double_t _time;
                Double_t _freq;
};

#endif
