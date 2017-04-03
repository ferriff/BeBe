#include "data_reader.h"

#include <cstdio>
#include <cstring>
#include <string>

#include "TFile.h"
#include "TObjString.h"
#include "TTree.h"

#include <libgen.h>


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
        TTree * t = _t.back();
        t->Branch("nsamples", &_nsamples, "nsamples/I");
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
        size_t cnt = 0;
        std::vector<std::string> tnames;
        for (auto t : _t) {
                tnames.push_back(t->GetName());
        }
        auto dir = gDirectory->GetPath();
        _fout->cd("MetaData");
        TTree * meta_data = new TTree("meta", "meta");
        meta_data->Branch("tree_order", &tnames);
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
        trigger_derivative(data, idx);
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
        std::vector<std::string> dnames;
        std::string header = "";
        char dname[128];
        _freq = 0;
        while ( (read = getline(&line, &len, fd)) != -1 ) {
                header += line;
                if (sscanf(line, "* Detecteur %s", &dname[0]) == 1) {
                        fprintf(stderr, "# New detector: `%s' --> ", dname);
                        if (_first_file) {
                                init_tree((std::string(dname) + "_heat").c_str());
                                init_tree((std::string(dname) + "_light").c_str());
                        }
                }
                if (sscanf(line, "Bolo.position = %x", &_detid) == 1) {
                        fprintf(stderr, "%d (detid, aka Bolo.position)\n", _detid);
                        // heat + light for each detector, binary data ordering
                        detids.push_back(_detid);        // heat
                        detids.push_back(_detid + 1000); // light
                }
                if (sscanf(line, "Echantillonage = %lf", &_freq) == 1) {
                        fprintf(stderr, "# Sampling frequency (kHz): %lf\n", _freq);
                        _freq *= 1000.;
                }
                if (strcmp(line, "* Donnees\n") == 0) {
                        fprintf(stderr, "# Header parsed.\n");
                        break;
                }
        }
        if (_ndetids && _ndetids != detids.size()) {
                fprintf(stderr, "[bb::data_reader::read_streamer_mode_file] Error: the header of `%s' is inconsistent with the previous ones:\n", input_file_name);
                fprintf(stderr, "[bb::data_reader::read_streamer_mode_file] ... found %lu detectors instead of %lu, aborting.\n", detids.size(), _ndetids);
                return;
        }

        // create association  detector name <-> detid
        if (_first_file) {
                for (size_t i = 0; i < dnames.size(); ++i) {
                        fprintf(stderr, "--> %s %d\n", dnames[i].c_str(), detids[i]);
                        //_detid_names.push_back(std::make_pair(dnames[i], detids[i]));
                        _detid_names[dnames[i]] = detids[i];
                }
        }
        // FIXME: add more checks on detid ordering when switching files...

        _ndetids = detids.size();
        _event_ids = (Long64_t *)calloc(_ndetids, sizeof(Long64_t));
        if (_first_file) {
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
        for (auto i = 0; i < _ndetids; ++i) trg[i] = 0;

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
                                        fprintf(stderr, "# Copying data: %lu --> %lu %lu    %d  %lu %lld\n", _cnt, idetid, idx, sample, _cnt / _ndetids, trg[idetid]);
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
                                        _t[idetid]->Fill();
                                        ++_event_ids[idetid];
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
