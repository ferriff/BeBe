#include "DataReader.h"

#include <cstdio>
#include <cstring>
#include <string>

#include "TFile.h"
#include "TObjString.h"
#include "TTree.h"

#include <libgen.h>


DataReader::DataReader(const char * output_file_name)
{
        _fout = TFile::Open(output_file_name, "recreate");
        if (!_fout) {
                fprintf(stderr, "[DataReader] problems in opening output file `%s', aborting\n", output_file_name);
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


DataReader::~DataReader()
{
        _t->Write(NULL, TObject::kWriteDelete);
        _fout->Close();
}


void DataReader::readStreamerModeFile(const char * input_file_name)
{
        fprintf(stdout, "# Reading file `%s'\n", input_file_name);
        FILE * fd = fopen(input_file_name, "r");
        if (!fd) {
                fprintf(stderr, "[DataReader::readStreamerModeFile] Cannot open file `%s'. Abort.", input_file_name);
                exit(3);
        }
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
                fprintf(stderr, "[DataReader::readFile] Error: the header of `%s' is inconsistent with the previous ones:\n", input_file_name);
                fprintf(stderr, "[DataReader::readFile] ... found %lu detectors instead of %lu, aborting.\n", detids.size(), _ndetids);
                exit(2);
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
        daqint * data[_ndetids];
        for (size_t i = 0; i < _ndetids; ++i) {
                data[i] = (daqint *)calloc(_nsamples, sizeof(daqint));
        }

        // read data and fill the ntuple
        daqint sample;
        int cnt = 0;
        while (fread(&sample, sizeof(daqint), 1, fd)) {
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
