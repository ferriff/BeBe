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
        // declare ntuple
        _t = new TTree("ntu", "ntu");
        _t->Branch("nsamples", &_nsamples, "nsamples/I");
        char tmp[64];
        sprintf(tmp, "raw_pulse[nsamples]/S");
        _br = _t->Branch("raw_pulse", (void *)0, tmp);
        _t->Branch("detid", &_detid, "detid/I");
        _t->Branch("time", &_time, "time/D");
        _t->Branch("freq", &_freq, "freq/D");
        _time = 0;
}


bb::data_reader::~data_reader()
{
        if (_t) _t->Write(NULL, TObject::kWriteDelete);
        if (_fout) _fout->Close();
}


void bb::data_reader::read_streamer_mode_file(const char * input_file_name)
{
        FILE * fd = fopen(input_file_name, "r");
        if (!fd) {
                fprintf(stderr, "[bb::data_reader::read_streamer_mode_file] Cannot open file `%s'. Abort.", input_file_name);
                return;
        }
        fprintf(stdout, "# Reading file `%s'\n", input_file_name);
	char * line = NULL;
	size_t len = 0;
	ssize_t read;

        // read and parse header
        std::vector<int> detids;
        std::string header = "";
        char dname[128];
        _freq = 0;
        while ( (read = getline(&line, &len, fd)) != -1 ) {
                header += line;
                if (sscanf(line, "* Detecteur %s", &dname[0]) == 1) {
                        fprintf(stderr, "# New detector: `%s' --> ", dname);
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
                fprintf(stderr, "[bb::data_reader::readFile] Error: the header of `%s' is inconsistent with the previous ones:\n", input_file_name);
                fprintf(stderr, "[bb::data_reader::readFile] ... found %lu detectors instead of %lu, aborting.\n", detids.size(), _ndetids);
                return;
        }
        _ndetids = detids.size();

        _fout->cd();
        TObjString os(header.c_str());
        char tmp[512];
        sprintf(tmp, "%s", input_file_name);
        std::string ttmp(basename(tmp));
        std::replace(ttmp.begin(), ttmp.end(), '/', '-');
        os.Write(("header_" + ttmp).c_str());

	// declare ntuple variables
        daqint_t * data[_ndetids];
        for (size_t i = 0; i < _ndetids; ++i) {
                data[i] = (daqint_t *)calloc(_nsamples, sizeof(daqint_t));
        }

        // read data and fill the ntuple
        daqint_t sample;
        int cnt = 0;
        while (fread(&sample, sizeof(daqint_t), 1, fd)) {
                if (cnt / (_ndetids * _nsamples) && cnt % (_ndetids * _nsamples) == 0) {
                        for (size_t i = 0; i < _ndetids; ++i) {
                                _br->SetAddress(data[i]);
                                _detid = detids[i];
                                _t->Fill();
                        }
                        _time += (Double_t)_nsamples / _freq;
                }
                data[cnt % _ndetids][(cnt / _ndetids) % _nsamples] = sample;
                ++cnt;
        }
        for (size_t i = 0; i < _ndetids; ++i) {
                free(data[i]);
        }
}


void bb::data_reader::read_trigger_mode_file(const char * heat_data_file, const char * light_data_file, const char * trigger_file)
{
        FILE * fd_heat = fopen(heat_data_file, "r");
        if (!fd_heat) {
                fprintf(stderr, "[bb::data_reader::read_trigger_mode_file] Cannot open file `%s'. Abort.", heat_data_file);
                exit(3);
        }
        FILE * fd_light = fopen(light_data_file, "r");
        if (!fd_light) {
                fprintf(stderr, "[bb::data_reader::read_trigger_mode_file] Cannot open file `%s'. Abort.", light_data_file);
                exit(3);
        }
        FILE * fd_trigger = fopen(trigger_file, "r");
        if (!fd_trigger) {
                fprintf(stderr, "[bb::data_reader::read_trigger_mode_file] Cannot open file `%s'. Abort.", trigger_file);
                exit(3);
        }
        fprintf(stderr, "[bb::data_reader::read_trigger_mode_file] Number of samples set to %lu. Make sure it is the correct one for the file read.", _nsamples);

        daqint_t * data;
        data = (daqint_t *)calloc(_nsamples, sizeof(daqint_t));
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
                        data[i] = (Int_t)sample;
                        //sample = bswap(sample);
                }
                _t->Fill();
                // read and fill light
                _detid += 1000;
                for (size_t i = 0; i < _nsamples; ++i) {
                        fread(&sample, sizeof(daqint_t), 1, fd_light);
                        data[i] = (Int_t)sample;
                        //sample = bswap(sample);
                }
                _t->Fill();
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
