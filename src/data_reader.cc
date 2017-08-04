#include "data_reader.h"

#include <cstdio>
#include <cstring>
#include <string>

#include "TFile.h"
#include "TObjString.h"
#include "TTree.h"

#include <libgen.h>


// trim from start (in place)
static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
            std::not1(std::ptr_fun<int, int>(std::isspace))));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
            std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}


void bb::data_reader::split(char * str, std::vector<std::string> & res)
{
        char * pch;
        pch = strtok (str," ");
        char tmp[128];
        while (pch != NULL) {
                sprintf(tmp, "%s", pch);
                res.push_back(tmp);
                pch = strtok(NULL, " ");
        }
}


bb::data_reader::data_reader(const char * output_file_name)
{
        _fout = TFile::Open(output_file_name, "recreate");
        if (!_fout || _fout->IsZombie()) {
                fprintf(stderr, "[bb::data_reader::data_reader] problems in opening output file `%s', aborting.\n", output_file_name);
                exit(-1);
        }
        _fout->mkdir("MetaData");
        _fout->mkdir("Data");
}


void bb::data_reader::init_tree(const char * name)
{
        // init ntuple
        _fout->cd("Data");
        _t.emplace_back(new TTree(name, name));
        _tree_order.push_back(name);
        _detids.push_back(_detid);
        TTree * t = _t.back();
        t->Branch("nsamples", &_nsamples, "nsamples/s");
        char tmp[64];
        sprintf(tmp, "raw_pulse[nsamples]/S");
        _br.emplace_back(t->Branch("raw_pulse", (void *)0, tmp));
        t->Branch("detid", &_detid, "detid/I");
        t->Branch("time", &_time, "time/D");
        t->Branch("freq", &_freq, "freq/D");
        _time = 0;
        _fout->cd();
}


void bb::data_reader::add_intra_tree_info()
{
        if (!_ndetids) {
                fprintf(stderr, "[bb::data_reader::add_intra_tree_info] Error: no detids found, make sure to call this function after the header has been parsed correctly.\n");
        }
        char tmp[64];
        sprintf(tmp, "event_ids[%lu]", _ndetids);
        for (auto t : _t) {
                t->Branch(tmp, _event_ids, (std::string(tmp) + "/L").c_str());
        }
}

void bb::data_reader::add_meta_data_info()
{
        auto dir = gDirectory->GetPath();
        _fout->cd("MetaData");
        TTree * meta_data = new TTree("meta", "meta");
        meta_data->Branch("tree_order", &_tree_order);
        meta_data->Branch("detid", &_detids);
        meta_data->Fill();
        meta_data->Write();
        delete meta_data;
        gDirectory->cd(dir);
}


void bb::data_reader::copy_data(size_t idx, daqint_t * d1, daqint_t * d2)
{
        std::copy(d1 + idx, d1 + _nsamples, d2);
        std::copy(d1, d1 + idx, d2 + (_nsamples - idx));
}


bb::data_reader::~data_reader()
{
        _fout->cd("Data");
        for (auto t : _t) {
                if (t) t->Write(NULL, TObject::kWriteDelete);
        }
        _fout->cd();
        if (_fout) _fout->Close();
        if (_event_ids) delete _event_ids;
}


bool bb::data_reader::trigger_over_threshold(daqint_t * data, size_t idx)
{
        return data[idx] > 1000;
}


bool bb::data_reader::trigger_over_threshold_with_baseline(daqint_t * data, size_t idx)
{
        size_t nped = _nsamples / 100;
        int start = (idx - nped);
        start = start > 0 ? start : _nsamples + start;
        float mean = 0;
        for (size_t i = 0; i < nped; ++i) mean += data[(start + i) % _nsamples];
        mean /= nped;
        //fprintf(stderr, "# mean: %f %d\n", mean, data[idx]);
        return data[idx] - mean > 1000;
}


bool bb::data_reader::trigger_derivative(daqint_t * data, size_t idx)
{
        return data[(idx + 25) % _nsamples] - data[idx] > 500;
}


bool bb::data_reader::trigger(daqint_t * data, size_t idx)
{
        //trigger_over_threshold_with_baseline(data, idx);
        return trigger_derivative(data, idx);
}


void bb::data_reader::read_streamer_mode_file(const char * input_file_name)
{
        FILE * fd = fopen(input_file_name, "r");
        if (!fd) {
                fprintf(stderr, "[bb::data_reader::read_streamer_mode_file] Cannot open file `%s'. Abort.\n", input_file_name);
                return;
        }
        fprintf(stdout, "# Reading file `%s'\n", input_file_name);
	char * line = NULL;
	size_t len = 0;
	ssize_t read;

        // read and parse header
        std::vector<int> detids;
        std::map<std::string, int> dnames;
        std::vector<std::vector<std::string> > chnames;
        std::vector<std::string> ch_unsplit_names;
        std::string header = "";
        char dname[128];
        _freq = 0;
        while ( (read = getline(&line, &len, fd)) != -1 ) {
                header += line;
                if (sscanf(line, "* Detecteur %s", &dname[0]) == 1) {
                        fprintf(stderr, "# New detector: `%s'", dname);
                }
                // attribute detids on the basis of the detector position
                // a unique identifier set by the DAQ
                if (sscanf(line, "Bolo.position = %x", &_detid) == 1) {
                        fprintf(stderr, "%d (detid, aka Bolo.position)\n", _detid);
                        // collect the detids from their position
                        dnames[dname] = _detid;
                }
                // read sampling frequency
                if (sscanf(line, "Echantillonage = %lf", &_freq) == 1) {
                        fprintf(stderr, "# Sampling frequency (kHz): %lf\n", _freq);
                        _freq *= 1000.;
                }
                // read the channel name
                // more than a channel can be associated to the same detector (heat, light, ionization, etc.)
                if (sscanf(line, "* Voie \"%[^\"]s\"", &dname[0]) == 1) {
                        fprintf(stderr, "# New channel: `%s'\n", dname);
                        ch_unsplit_names.push_back(dname);
                        std::vector<std::string> res;
                        split(dname, res);
                        chnames.push_back(res);
                }
                if (strcmp(line, "* Donnees\n") == 0) {
                        fprintf(stderr, "# Header parsed.\n");
                        break;
                }
        }
        // for each channel, look for the associated detector and build a unique identifier
        // made from the detid and the channel sequential number
        size_t cnt_tmp(0);
        std::map<int, int> tmpm;
        for (auto & ch : chnames) {
                for (auto & n : ch) {
                        auto it = dnames.find(n);
                        if (it != dnames.end()) {
                                auto nit = ch_unsplit_names[cnt_tmp].find(n);
                                if (nit != std::string::npos) {
                                        ch_unsplit_names[cnt_tmp].erase(nit, n.size());
                                        trim(ch_unsplit_names[cnt_tmp]);
                                }
                                auto mit = tmpm.find(it->second);
                                if (mit == tmpm.end()) tmpm[it->second] = 0;
                                else tmpm[it->second] += 1;
                                _detid = it->second * 100 + mit->second;
                                fprintf(stderr, "Found detid %d for channel `%s' of detector `%s', attributing detid %d\n", it->second, ch_unsplit_names[cnt_tmp].c_str(), it->first.c_str(), _detid);
                                // if the file is the first read, initialize the trees
                                if (_first_file) {
                                        init_tree((it->first + "_" + ch_unsplit_names[cnt_tmp]).c_str());
                                }
                                detids.push_back(_detid);
                        }
                }
                ++cnt_tmp;
        }

        if (_ndetids && _ndetids != detids.size()) {
                fprintf(stderr, "[bb::data_reader::read_streamer_mode_file] Error: the header of `%s' is inconsistent with the previous ones:\n", input_file_name);
                fprintf(stderr, "[bb::data_reader::read_streamer_mode_file] ... found %lu detectors instead of %lu, aborting.\n", detids.size(), _ndetids);
                return;
        }

        _ndetids = detids.size();
        if (_first_file) {
                _event_ids = (Long64_t *)calloc(_ndetids, sizeof(Long64_t));
                for (size_t i = 0; i < _ndetids; ++i) _event_ids[i] = -1;
                add_intra_tree_info();
                add_meta_data_info();
        }

        // write header information in the output file
        _fout->cd("MetaData");
        TObjString os(header.c_str());
        char tmp[512];
        sprintf(tmp, "%s", input_file_name);
        std::string ttmp(basename(tmp));
        std::replace(ttmp.begin(), ttmp.end(), '/', '-');
        os.Write(("header_" + ttmp).c_str());
        _fout->cd();

        // write the ordering of the trees, to reconstruct the association
        //    position in the vector of indices <-> tree
        _fout->cd("MetaData");
        _fout->cd();

	// declare ntuple variables
        daqint_t * data[_ndetids];
        for (size_t i = 0; i < _ndetids; ++i) {
                data[i] = (daqint_t *)calloc(_nsamples, sizeof(daqint_t));
        }
        daqint_t cdata[_nsamples];
        daqint_t sample[_ndetids];

        // read data and fill the ntuple
        //  1)  fill an array of data
        //      det1 sample1, s2, s3, ..., s_nsamples
        //      ...                               ...
        //      detN sample1, s2, s3, ..., s_nsamples
        //  2)  if the sample S for the detector D is above a trigger threshold
        //      then
        //        . write the event in the corresponding tree
        //        . write the tree indices of the previous events of all the
        //          other trees (to ease navigation)
        //size_t ntot = _ndetids * _nsamples;

        // to hold the last triggered sample
        Long64_t trg[_ndetids];
        for (size_t i = 0; i < _ndetids; ++i) trg[i] = 0;

        // read an array of ndetids_
        while (fread(&sample, _ndetids * sizeof(daqint_t), 1, fd)) {
                // start looking for triggers when the array is fully filled
                // from 1/3 of the array, i.e. the first 1/3 is kept for
                // the pre-samples
                /////////////fprintf(stderr, "%d  ", nb);
                /////////////for (size_t i = 0; i < _ndetids; ++i) fprintf(stderr, " %10d", sample[i]);
                /////////////fprintf(stderr, "\n");
                for (size_t idetid = 0; idetid < _ndetids; ++idetid) {
                        if (_cnt > _nsamples) {
                                // if sample above trigger threshold
                                // and outside the trigger window, save the event
                                // in the tree of the corresponding detid
                                //size_t idetid = _cnt % _ndetids;
                                // divide _nsamples in 3 equal parts:
                                // 1. pre-samples (no triggers searched for)
                                // 2. one and only one trigger allowed
                                // 3. post-samples (new triggers allowed)
                                size_t pre_samples = _nsamples / 3;
                                size_t core_samples = _nsamples / 3;
                                // start looking for triggers from idx
                                size_t idx = (_cnt + pre_samples) % _nsamples;
                                if (trigger(data[idetid], idx) && (_cnt - trg[idetid]) > core_samples) {
                                        ++_event_ids[idetid];
                                        trg[idetid] = _cnt;
                                        int start = (idx - pre_samples);
                                        bzero(cdata, sizeof(cdata));
                                        copy_data(start > 0 ? start : _nsamples + start, data[idetid], cdata);
                                        //if (idetid == 6) {
                                        //        for (size_t i = 0; i < _nsamples; ++i)
                                        //        {
                                        //                fprintf(stderr, "%lu %d %d\n", i, data[idetid][i], cdata[i]);
                                        //        }
                                        //        fprintf(stderr, "\n\n");
                                        //}
                                        _br[idetid]->SetAddress(cdata);
                                        _detid = detids[idetid];
                                        fprintf(stderr, "# Copying data: %10lu --> %d %2lu %4.4f %5lu    %hn  %8lu %8lld %8lld", _cnt, detids[idetid], idetid, _time, idx, sample, _cnt / _ndetids, trg[idetid], _event_ids[idetid]);
                                        fprintf(stderr, "  ");
                                        for (size_t i = 0; i < _ndetids; ++i) {
                                                fprintf(stderr, " %lld", _event_ids[i]);
                                        }
                                        fprintf(stderr, "\n");
                                        _t[idetid]->Fill();
                                }
                                //// keep track of the time
                                //_time += (Double_t)_nsamples / _freq;
                        }
                        // keep track of the time
                        _time = (Double_t)(_cnt) / _freq;
                        // fill the data array
                        data[idetid][_cnt % _nsamples] = sample[idetid];
                }
                ++_cnt;
        }
        for (size_t i = 0; i < _ndetids; ++i) {
                free(data[i]);
        }
        _first_file = false;
}


void bb::data_reader::read_trigger_mode_file(const char * heat_data_file, const char * light_data_file, const char * trigger_file)
{
        FILE * fd_heat = fopen(heat_data_file, "r");
        if (!fd_heat) {
                fprintf(stderr, "[bb::data_reader::read_trigger_mode_file] Cannot open file `%s'. Abort.\n", heat_data_file);
                exit(3);
        }
        FILE * fd_light = fopen(light_data_file, "r");
        if (!fd_light) {
                fprintf(stderr, "[bb::data_reader::read_trigger_mode_file] Cannot open file `%s'. Abort.\n", light_data_file);
                exit(3);
        }
        FILE * fd_trigger = fopen(trigger_file, "r");
        if (!fd_trigger) {
                fprintf(stderr, "[bb::data_reader::read_trigger_mode_file] Cannot open file `%s'. Abort.\n", trigger_file);
                exit(3);
        }
        fprintf(stderr, "[bb::data_reader::read_trigger_mode_file] Number of samples set to %lu. Make sure it is the correct one for the file read.\n", _nsamples);

        init_tree("ntu");

        daqint_t * data;
        data = (daqint_t *)calloc(_nsamples, sizeof(daqint_t));
        _br[0]->SetAddress(data);
        daqint_t sample;

        char * line = NULL;
        size_t len = 0;
        ssize_t read;
        size_t is = 0;

        while ( (read = getline(&line, &len, fd_trigger)) != -1 ) {
                if (line[0] == '#') continue;
                sscanf(line, "%*d %*f %*f %d %lf", &_detid, &_time);
                // read and fill heat
                for (size_t i = 0; i < _nsamples; ++i) {
                        fread(&sample, sizeof(daqint_t), 1, fd_heat);
                        //sample = bswap(sample);
                        data[i] = sample;
                }
                _t[0]->Fill();
                // read and fill light
                _detid += 1000;
                for (size_t i = 0; i < _nsamples; ++i) {
                        fread(&sample, sizeof(daqint_t), 1, fd_light);
                        //sample = bswap(sample);
                        data[i] = sample;
                }
                _t[0]->Fill();
                ++is;
        }
        if (data) free(data);
        if (fread(&sample, sizeof(daqint_t), 1, fd_heat) != 0 || fread(&sample, sizeof(daqint_t), 1, fd_light) != 0) {
                fprintf(stderr, "[bb::data_reader::read_trigger_mode_file] ERROR: data"
                        " files seem to contain more events than"
                        " trigger file, abort\n");
                exit(1);
        }
}
