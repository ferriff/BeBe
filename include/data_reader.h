#ifndef DATA_READER_H
#define DATA_READER_H
/*
 *        Class: data_reader
 *  Description: handling of SAMBA binarary file in streamer mode
 *       Author: Federico Ferri, CEA/Saclay
 */
#include "types.h"

#include "TFile.h"
#include "TObjString.h"
#include "TTree.h"

namespace bb {

        class data_reader
        {
                public:
                        data_reader(const char * output_file_name);

                        ~data_reader();

                        size_t n_samples() const { return _nsamples; }
                        void set_n_samples(size_t ns) { _nsamples = ns; }

                        void read_streamer_mode_file(const char * input_file_name);
                        void read_trigger_mode_file(const char * heat_data_file, const char * light_data_file, const char * trigger_file);

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

}

#endif
