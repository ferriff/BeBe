#include <byteswap.h>
#include <float.h>
#include <stdio.h>
#include <stdint.h>

#include "TFile.h"
#include "TTree.h"

#define FIFO_LENGTH          500
#define TIME_WINDOW         2000
#define PRE_TRIGGER          400 // must be lower than FIFO_LENGTH
#define TRIGGER_THRESHOLD   2000
#define PED_SAMPLES           50 // number of samples on which the pedestal is computed
#define PED               -80184
#define ADCTOEV            21. / (1<<18 - 1) / 165.e-9 / 1000. // dynamics (V) / bits / (165 Volts / keV)

typedef struct {
        size_t idx;
        int32_t data[FIFO_LENGTH];
} fifo;

void fifo_reset(fifo * f)
{
        f->idx = 0;
        for (size_t i = 0; i < FIFO_LENGTH; ++i) f->data[i] = 0;
}

void fifo_in(fifo * f, int32_t value)
{
        f->data[f->idx++ % FIFO_LENGTH] = value;
}

int32_t fifo_read(const fifo * f, size_t pos)
{
        return f->data[(f->idx + pos) % FIFO_LENGTH];
}

int32_t signal_max(const int32_t data[], size_t size, size_t * idx)
{
        int32_t m = INT32_MIN;
        size_t im = 0;
        for (size_t i = 0; i < size; ++i) {
                if (data[i] > m) {
                        m = data[i];
                        im = i;
                }
        }
        if (idx) *idx = im;
        return m;
}

float signal_average(const int32_t data[], size_t size)
{
        int32_t m = 0;
        for (size_t i = 0; i < size; ++i) {
                m += data[i];
        }
        return (float)m / (float)size;
}

int32_t signal_quality(const int32_t data[])
{
        for (size_t i = 300; i < TIME_WINDOW; ++i) {
                if (data[i] - data[i - PED_SAMPLES] > 0) return 0;
        }
        return 1;
}

float signal_convert(const int32_t data[], float fata[])
{
        float ped = signal_average(data, PED_SAMPLES);
        for (size_t i = 0; i < TIME_WINDOW; ++i) fata[i] = data[i] - ped;
        return ped;
}

float fignal_max(const float data[], size_t size, size_t * idx)
{
        float m = FLT_MIN;
        size_t im = 0;
        for (size_t i = 0; i < size; ++i) {
                if (data[i] > m) {
                        m = data[i];
                        im = i;
                }
        }
        if (idx) *idx = im;
        return m;
}

size_t fignal_decay_time(const float data[], size_t imax, float fraction)
{
        float m = data[imax];
        for (size_t i = imax; i < TIME_WINDOW; ++i) {
                if (data[i] / m < fraction) return i - imax;
        }
        return 0;
}

size_t fignal_rise_time(const float data[], size_t imax, float fraction)
{
        float m = data[imax];
        for (size_t i = imax; i >= 0; --i) {
                if (data[i] / m < fraction) return imax - i;
        }
        return 0;
}

size_t pileup(const float data[], size_t imax)
{
        for (size_t i = imax; i < TIME_WINDOW - 50; ++i) {
                if (data[i + 10] - data[i] > 50) {
                        //for (size_t j = 0; j < TIME_WINDOW; ++j) {
                        //        fprintf(stderr, "%d %f\n", j, data[j]);
                        //}
                        return i + 10;
                }
        }
        return 0;
}

void read_data(FILE * fd)
{
        if (!fd) return;
        fifo f;
        fifo_reset(&f);
        char trigger = 0;
        size_t nsamples = TIME_WINDOW;
        int32_t * data = (int32_t *)calloc(nsamples, sizeof(int32_t));
        float   * fata = (float *)calloc(nsamples, sizeof(float));
        int32_t sample;
        int32_t detid = 0, pu = 0;
        size_t idx = 0;
        size_t tdaq = 0, t_max_daq = 0, t_start_daq = 0;
        size_t im = 0, im_tmp = 0, im_prev = 0;
        float m = 0., p = 0.;
        size_t decay_time = 0, rise_time = 0;
        TFile * fout = TFile::Open("signals.root", "recreate");
        TTree * t = new TTree("ntu", "ntu");
        char tmp[64];
        //sprintf(tmp, "raw_signal[%d]/I", TIME_WINDOW);
        t->Branch("nsamples", &nsamples, "nsamples/i");
        sprintf(tmp, "raw_signal[nsamples]/I");
        TBranch * br = t->Branch("raw_signal", data, tmp);
        t->Branch("t_max", &im, "t_max/i");
        t->Branch("t_max_daq", &t_max_daq, "t_max_daq/i");
        t->Branch("t_start_daq", &t_start_daq, "t_start_daq/i");
        t->Branch("amplitude", &m, "amplitude/F");
        t->Branch("baseline", &p, "baseline/F");
        t->Branch("detid", &detid, "detid/I");
        t->Branch("pu", &pu, "pu/I");
        while (fread(&sample, sizeof(int32_t), 1, fd)) {
                sample = bswap_32(sample);
                sample -= PED;
                fifo_in(&f, sample);
                if (sample > TRIGGER_THRESHOLD && !trigger) {
                        for (size_t i = 0; i < PRE_TRIGGER; ++i) data[idx++] = fifo_read(&f, -PRE_TRIGGER + i);
                        trigger = 1;
                        --idx;
                }
                if (trigger) {
                        data[idx++] = sample;
                        if (idx == TIME_WINDOW) {
                                trigger = 0;
                                //if (!signal_quality(data)) continue;
                                p = signal_convert(data, fata);
                                /* integer analysis
                                size_t im = 0;
                                int32_t ped = signal_average(data, PED_SAMPLES);
                                int32_t m = signal_max(data, 300, &im) - ped;
                                if (m < 2000) continue;
                                for (int i = 0; i < TIME_WINDOW; ++i) fprintf(stdout, "%d %d %lu %d\n", i, data[i] - ped, im, m);
                                */
                                im = 0;
                                m = fignal_max(fata, PRE_TRIGGER + 100, &im);
                                decay_time = fignal_decay_time(fata, im, 0.3);
                                rise_time  = fignal_rise_time(fata, im, 0.3);
                                //fprintf(stderr, "%lu %d %lu\n", tdaq, TIME_WINDOW, im);
                                //getchar();
                                t_start_daq = tdaq > TIME_WINDOW ? tdaq - TIME_WINDOW : 0;
                                t_max_daq = tdaq > TIME_WINDOW ? tdaq - TIME_WINDOW + im : 0;
                                nsamples = TIME_WINDOW;
                                // detect pileup
                                im_tmp = im;
                                im_prev = im;
                                pu = 0;
                                while ((im_tmp = pileup(fata, im_tmp)) != 0) {
                                        nsamples += im_tmp - im_prev;
                                        data = (int32_t *)realloc(data, nsamples * sizeof(int32_t));
                                        br->SetAddress(data);
                                        fprintf(stderr, "nsamples = %lu  im_tmp = %lu  im_prev = %lu (im = %lu)  sizeof(data) = %lu\n", nsamples, im_tmp, im_prev, im, sizeof(data) / sizeof(int32_t));
                                        im_prev = im_tmp;
                                        while (fread(&sample, sizeof(int32_t), 1, fd) && idx < nsamples) {
                                                sample = bswap_32(sample);
                                                sample -= PED;
                                                data[idx++] = sample;
                                                //fprintf(stderr, "--> %lu  %d %d\n", idx, sample, data[idx - 1]);
                                        }
                                        pu += 1;
                                }
                                idx = 0;
                                t->Fill();
                                //if (m < 2000.) continue;
                                //fprintf(stdout, "# %d %lu %lu %f %lu %lu\n", TIME_WINDOW, tdaq, im, m, rise_time, decay_time);
                                //for (int i = 0; i < TIME_WINDOW; ++i) fprintf(stdout, "%d %f %lu %f %lu %lu\n", i, fata[i], im, m, rise_time, decay_time);
                                //fprintf(stdout, "\n\n");
                        }
                }
                ++tdaq;
        }
        free(data);
        free(fata);
        t->Write();
        fout->Close();
}

int main()
{
        FILE * fd = fopen("Raw_501311.BIN17", "r");
        read_data(fd);
        return 0;
}
