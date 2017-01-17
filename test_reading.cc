#include <cstdio>
#include <cstring>
#include <string>

#include <vector>

#include "TFile.h"
#include "TObjString.h"
#include "TTree.h"

typedef int16_t daqint;
#define bswap(sample) { bswap_16(sample) }


int main()
{
        char filename [] = "../../data/betabeta/ra13e000_DATA_FOR_TEST/Raw_Modane_stream/ra13e000_000";
        FILE * fd = fopen(filename, "r");
	char * line = NULL;
	size_t len = 0;
	ssize_t read;

        // declare output file
        TFile * fout = TFile::Open("pappa.root", "recreate");
        fout->cd();

        // read and parse header
        int32_t detid;
        std::vector<int> detids;
        std::string header = "";
        while ( (read = getline(&line, &len, fd)) != -1 ) {
                header += line;
                if (sscanf(line, "Bolo.position = %x", &detid) == 1) {
                        fprintf(stderr, "# Found new detector in position: %d\n", detid);
                        detids.push_back(detid);        // heat
                        detids.push_back(detid + 1000); // light
                }
                if (strcmp(line, "* Donnees\n") == 0) {
                        fprintf(stderr, "# Header parsed.\n");
                        break;
                }
        }
        // double the channels (heat + light for each detector)
        //for (auto v : detids) detids.push_back(v + 1000);
        size_t ndetids = detids.size();

        TObjString os(header.c_str());
        os.Write("header");

	// declare ntuple variables
        size_t nsamples = 1000;
        daqint * data[ndetids];
        for (size_t i = 0; i < ndetids; ++i) {
                data[i] = (daqint *)calloc(nsamples, sizeof(daqint));
        }
        //float   * fata = (float *)calloc(nsamples, sizeof(float));
        // declare ntuple
        TTree * t = new TTree("ntu", "ntu");
        t->Branch("nsamples", &nsamples, "nsamples/I");
        char tmp[64];
        sprintf(tmp, "raw_pulse[nsamples]/S");
        TBranch * br = t->Branch("raw_pulse", data[0], tmp);
        t->Branch("detid", &detid, "detid/I");
        //t->Branch("time", &ts, "time/D");

        // read data and fill ntuple
        daqint sample;
        int cnt = 0;
        while (fread(&sample, sizeof(daqint), 1, fd)) {
                if (cnt / (ndetids * nsamples) && cnt % (ndetids * nsamples) == 0) {
                        for (size_t i = 0; i < ndetids; ++i) {
                                //fprintf(stdout, "filling %d (%d) [%lu]\n", detids[i], cnt, nsamples);
                                br->SetAddress(data[i]);
                                detid = detids[i];
                                t->Fill();
                        }
                }
                //detid = detids[(cnt / 2) % ndetids] + 1000 * (cnt % 2);
                data[cnt % ndetids][(cnt / ndetids) % nsamples] = sample; // all heats then all lights
                //fprintf(stderr, "detid = %d  sample = %d\n", , sample);
                //fprintf(stderr, "detid = %d  sample = %d\n", detids[cnt % ndetids] + 1000 * ((cnt / ndetids) % 2), sample);
                //if (cnt > 500000) break;
                ++cnt;
        }
        t->Write(NULL, TObject::kWriteDelete);
        fout->Close();
        return 0;
}
