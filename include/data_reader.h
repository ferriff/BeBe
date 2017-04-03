#ifndef DATA_READER_H
#define DATA_READER_H
/*
 *        Class: data_reader
 *  Description: handling of SAMBA binarary file in trigger and streamer mode
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

                        void init_tree(const char * name);
                        void add_intra_tree_info();
                        void add_meta_data_info();
                        void copy_data(size_t idx, daqint_t * d1, daqint_t * d2);

                        bool trigger(daqint_t * data, size_t idx);
                        bool trigger_derivative(daqint_t * data, size_t idx);
                        bool trigger_over_threshold(daqint_t * data, size_t idx);
                        bool trigger_over_threshold_with_baseline(daqint_t * data, size_t idx);

                        std::vector<TTree *> _t;
                        std::vector<TBranch *> _br;
                        TFile * _fout;
                        FILE *  _fd;
                        size_t _nsamples = 3000;
                        size_t _ndetids = 0;
                        size_t _cnt = 0;
                        int32_t _detid;
                        std::map<std::string, int32_t> _detid_names;
                        Double_t _time;
                        Double_t _freq;
                        Long64_t * _event_ids;
                        bool _first_file = true;
        };

}

#endif
