#include <cstdio>

#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLeaf.h"
#include "TProfile.h"
#include "TTree.h"

#define TIME_WINDOW 2000
#define PRE_TRIGGER  400
#define PED_SAMPLES  100

typedef struct {
        int32_t detid;
        float   ampl;
        float   bl;
        UInt_t  t_max;
        UInt_t  t_max_daq;
        UInt_t  t_start_daq;
        UInt_t  pu;
        UInt_t  nsamples;
        int32_t * data;
} Event;

Event  _e;
size_t _ipulse;

void set_branches(TTree * t)
{
        t->SetBranchAddress("detid", &_e.detid);
        // handle variable size array
        // allocate maximum dimension (cf. TTreePlayer::MakeClass code)
        TLeaf * l = t->GetLeaf("nsamples", "nsamples");
        _e.data = (int32_t *)calloc(l->GetMaximum(), sizeof(int32_t));
        t->SetBranchAddress("raw_signal", _e.data);
        t->SetBranchAddress("amplitude", &_e.ampl);
        t->SetBranchAddress("baseline", &_e.bl);
        t->SetBranchAddress("t_max", &_e.t_max);
        t->SetBranchAddress("t_max_daq", &_e.t_max_daq);
        t->SetBranchAddress("t_start_daq", &_e.t_start_daq);
        t->SetBranchAddress("pu", &_e.pu);
        t->SetBranchAddress("nsamples", &_e.nsamples);
        //return br;
}


float sample_average(const int32_t data[], size_t start, size_t size)
{
        int32_t m = 0;
        for (size_t i = start; i < size; ++i) {
                m += data[i];
        }
        return (float)m / (float)size;
}


float sample_rms(const int32_t data[], size_t start, size_t size)
{
        float m = 0, mm = 0;
        for (size_t i = start; i < size; ++i) {
                m  += data[i];
                mm += data[i] * data[i];
        }
        m /= (float)(size - start);
        return sqrt(mm / (float)(size - start) - m * m);
}


void signal_convert(const int32_t data[], float fata[])
{
        float ped = sample_average(data, 0, PED_SAMPLES);
        for (size_t i = 0; i < TIME_WINDOW; ++i) fata[i] = data[i] - ped;
}


size_t signal_decay_time(const float data[], size_t imax, float fraction)
{
        float m = data[imax];
        for (size_t i = imax; i < TIME_WINDOW; ++i) {
                if (data[i] / m < fraction) return i - imax;
        }
        return 0;
}


size_t signal_rise_time(const float data[], size_t imax, float fraction)
{
        float m = data[imax];
        for (size_t i = imax; i >= 0; --i) {
                if (data[i] / m < fraction) return imax - i;
        }
        return 0;
}


float signal_rise_time_interpolated(const float data[], size_t imax, float fraction, float amplitude, float tmax)
{
        float m = data[imax];
        if (amplitude) m = amplitude;
        for (size_t i = imax; i >= 0; --i) {
                if (data[i] / m < fraction) {
                        //fprintf(stderr, "===> %f %f\n", imax - i + (m * fraction - data[i]) / (data[i + 1] - data[i]) * 1., imax - (i + 1 + (data[i + 1] - m * fraction) / (data[i + 1] - data[i]) * 1.));
                        //return imax - i + (m * fraction - data[i]) / (data[i + 1] - data[i]) * 1.;
                        return tmax - (i + 1 + (data[i + 1] - m * fraction) / (data[i + 1] - data[i]) * 1.);
                }
        }
        return -1;
}


float signal_decay_time_interpolated(const float data[], size_t imax, float fraction, float amplitude, float tmax)
{
        float m = data[imax];
        if (amplitude) m = amplitude;
        for (size_t i = imax; i < TIME_WINDOW; ++i) {
                if (data[i] / m < fraction) {
                        return i - (m * fraction - data[i]) / (data[i - 1] - data[i]) - tmax;
                }
        }
        return -1;
}


float signal_max(const float data[], size_t size, size_t * idx)
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


float signal_max_fitted(const float data[], size_t size, float * xM)
{
        float m = FLT_MIN;
        size_t im = 0;
        for (size_t i = 0; i < size; ++i) {
                if (data[i] > m) {
                        m = data[i];
                        im = i;
                }
        }
        char tmp[16];
        sprintf(tmp, "pulse_%06lu", _ipulse);
        TGraph * g = new TGraph();
        g->SetNameTitle(tmp, tmp);
        int npoints = 10;
        for (int i = 0; i < npoints; ++i) g->SetPoint(i, (float)im - npoints / 2 + i, (float)data[im - npoints / 2 + i]);
        g->Fit("pol2", "Q");
        TF1 * f = g->GetFunction("pol2");
        float p0 = f->GetParameter(0);
        float p1 = f->GetParameter(1);
        float p2 = f->GetParameter(2);
        if (xM) *xM = -0.5 * p1 / p2;
        delete g;
        return p0 - 0.25 * p1 * p1 / p2;
}


int signal_flag(float ampl, float trise, float tdecay, float pedrms)
{
//triplet = "$13 < 10 && $8 > 30 && $8 < 34 && $10 / 1000. > 3 && $10 / 1000. < 10"
//top1    = "$13 < 10 && $8 > 43 && $8 < 44 && $10 / 1000 < 30"
//top2    = "$13 < -10 && $8 > 23.25 && $8 < 24 && $10 / 1000 > 4"
//top3    = "$13 < 10 && $8 > 41.1 + 0.02 * $10 / 1000. && $8 < 42.5 + 0.01 * $10 / 1000. && $10 / 1000 < 60"
//top4    = "$13 < 10 && $8 > 40.25 + 0.0075 * $10 / 1000.  && $8 < 40.85 + 0.015 * $10 / 1000. && $10 / 1000 > 15"
//top5    = "$13 < 10 && $8 > 38.75 + 0.01 * $10 / 1000. && $8 < 39.5 + .01 * $10 / 1000. && $10 / 1000 > 10"
//top6    = "$13 < 10 && $8 > 38 && $8 < 39 && $10 / 1000 > 90"
//heat    = "$13 < 10 && $8 > 36  && $8 < 39 && $10 / 1000 > 40 && $10 / 1000. < 50"
        if (pedrms > 10) return -1;
        if (0);
        else if (trise > 30 && trise < 34 && ampl > 3e+3 && ampl < 1e+4) return 8;
        else if (trise > 43 && trise < 44 && ampl < 3e+4)                return 2;
        else if (trise > 41.1 + 0.02 * ampl / 1000. && trise < 42.5 + 0.01 * ampl / 1000. && ampl < 6e+4) return 3;
        else if (trise > 40.25 + 0.0075 * ampl / 1000.  && trise < 40.85 + 0.015 * ampl / 1000. && ampl > 15e+3) return 4;
        else if (trise > 38.75 + 0.01 * ampl / 1000. && trise < 39.5 + .01 * ampl / 1000. && ampl > 1e+4) return 5;
        else if (trise > 38 && trise < 39 && ampl > 9e+4) return 9;
        // heater
        else if (trise > 36  && trise < 39 && ampl > 4e+4 && ampl < 5e+4) return 7;
        return 0;
}


int main()
{
        _ipulse = 0;
        TFile * fout = TFile::Open("histos.root", "recreate");
        TH1F * h_ampl_raw               = new TH1F("h_raw_ampl", "h_raw_ampl", 2<<17, 0., 2<<17);
        TH1F     * h_ampl                   = (TH1F*)h_ampl_raw->Clone("h_ampl");
        TH1F     * h_ampl_drift_corr        = (TH1F*)h_ampl_raw->Clone("h_ampl_drift_corr");
        TH1F     * h_ampl_heater            = (TH1F*)h_ampl_raw->Clone("h_ampl_heater");
        TH1F     * h_ampl_heater_drift_corr = (TH1F*)h_ampl_raw->Clone("h_ampl_heater_drift_corr");
        TH1F     * h_ampl_keV               = new TH1F("h_ampl_keV", "h_ampl_keV", 4000, 0., 10000.);
        TH1F     * h_ped                    = new TH1F("h_ped", "h_ped", 2<<13, 0., 2<<17);
        TH1F     * h_ped_rms                = new TH1F("h_ped_rms", "h_ped_rms", 100, 0., 100.);
        TH1F     * h_dt                     = new TH1F("h_dt", "h_dt", 5000, 0., 5000000.);
        TH1F     * h_time_s1                = new TH1F("h_time_s1", "h_time_s1", 20000, 0., 50.);
        TProfile * p_average_signal         = new TProfile("p_average_signal", "p_average_signal", 5. * TIME_WINDOW, 0., TIME_WINDOW);
        TGraph * g_ampl_vs_t            = new TGraph();
        g_ampl_vs_t->SetNameTitle("g_ampl_vs_t", "g_ampl_vs_t");
        gDirectory->Add(g_ampl_vs_t);
        TGraph * g_ampl2_vs_t = new TGraph();
        g_ampl2_vs_t->SetNameTitle("g_ampl2_vs_t", "g_ampl2_vs_t");
        gDirectory->Add(g_ampl2_vs_t);
        TGraph * g_ampl_all_vs_t        = new TGraph();
        g_ampl_all_vs_t->SetNameTitle("g_ampl_all_vs_t", "g_ampl_all_vs_t");
        gDirectory->Add(g_ampl_all_vs_t);
        TFile * fin = TFile::Open("signals.root");
        TTree * t = (TTree*)fin->Get("ntu");
        fout->cd();
        set_branches(t);
        Long64_t nentries = t->GetEntries();
        Long64_t gcnt = 0, gcnt2 = 0, gcnta = 0;
        UInt_t otmaxdaq = 1;
        float fata[TIME_WINDOW]; // FIXME
        FILE * fp = fopen("pulses.dat", "w");
        for (Long64_t ien = 0; ien < nentries; ++ien) {
                t->GetEntry(ien);
                h_ampl_raw->Fill(_e.ampl);
                // signal analysis
                signal_convert(_e.data, fata);
                float ped = sample_average(_e.data, 0, 50);
                float ped_rms = sample_rms(_e.data, 0, 50);
                h_ped->Fill(ped);
                h_ped_rms->Fill(ped_rms);
                size_t iM = 0;
                float M = signal_max(fata, PRE_TRIGGER + 100, &iM);
                float fiM = -1;
                float fM = signal_max_fitted(fata, PRE_TRIGGER + 100, &fiM);
                float ft_daq = _e.t_max_daq - (_e.t_max - iM);
                float trise = signal_rise_time_interpolated(fata, iM, 0.05, fM, fiM) - signal_rise_time_interpolated(fata, iM, 0.95, fM, fiM);
                float tdecay = signal_decay_time_interpolated(fata, iM, 0.20, fM, fiM) - signal_decay_time_interpolated(fata, iM, 0.95, fM, fiM);
                fprintf(stdout, "%f %f %lu %u %u %lu %lu %f %f %f %f %f %f\n",
                        _e.ampl, M, iM, _e.t_max, _e.t_max_daq, 
                        signal_rise_time(fata, iM, 0.05), signal_decay_time(fata, iM, 0.20),
                        trise,
                        tdecay, 
                        fM, fiM, ped, ped_rms);
                h_ampl->Fill(fM);
                float adc2keV = 4783. / 101279.;
                h_ampl_keV->Fill(fM * adc2keV);
                h_ampl_drift_corr->Fill(fM - 1.788 * ft_daq / 2000. / 3600.);
                //if (otmaxdaq && _e.ampl > 101100 && _e.ampl < 101500) {
                if (otmaxdaq && _e.ampl > 46400 && _e.ampl < 46500) {
                        h_dt->Fill(_e.t_max_daq - otmaxdaq);
                        h_ampl_heater->Fill(fM);
                        h_ampl_heater_drift_corr->Fill(fM - 1.788 * ft_daq / 2000. / 3600.);
                        otmaxdaq = _e.t_max_daq;
                }
                if (_e.ampl > 46400 && _e.ampl < 46500) {
                        g_ampl_vs_t->SetPoint(gcnt++, ft_daq / 2000. / 3600., fM);
                } else if (_e.ampl > 100.8e+3 && _e.ampl < 101.6e+3) {
                        g_ampl2_vs_t->SetPoint(gcnt2++, ft_daq / 2000. / 3600., fM);
                }
                g_ampl_all_vs_t->SetPoint(gcnta++, (float)_e.t_start_daq / 2000. / 3600., fM);
                if (fM > 1e5) h_time_s1->Fill(trise - 0.01 * fM / 1000.);
                //else if (_e.ampl > 101.2e+3 && _e.ampl < 101.4e+3) g_ampl2_vs_t->SetPoint(gcnt2++, ft_daq / 2000. / 3600., fM);
                int flag = signal_flag(fM, trise, tdecay, ped_rms);
                if (fM > 1e4) {
                        fprintf(fp, "#flag: %d\n", flag);
                        for (size_t i = 0; i < TIME_WINDOW; ++i) {
                                fprintf(fp, "%lu %f %f %f\n", i, 450 - fiM, fata[i], fM);
                        }
                        fprintf(fp, "\n\n");
                }
                if (_e.pu == 0 && fM > 101e+3 && fM < 101.75e+3 && ien > 750 && ien < 1000) {
                        for (size_t is = 0; is < TIME_WINDOW; ++is) {
                                //fprintf(stderr, "--> %lu %lu %f %f %f\n", is, iM, fiM, fata[is], fM );
                                //if (is % 100 == 0) getchar();
                                p_average_signal->Fill(is - (fiM - 450), fata[is] / fM);
                        }
                }
                ++_ipulse;
        }
        if (_e.data) free(_e.data);
        h_ampl_drift_corr->Rebin(16);
        h_ampl_heater->Rebin(6);
        h_ampl_heater_drift_corr->Rebin(6);
        TH1F * h_ampl_rebin = (TH1F*)h_ampl->Clone("h_ampl_rebin");
        h_ampl_rebin->Rebin(12);
        h_ampl_rebin->GetXaxis()->SetRangeUser(100e+3, 102e+3);
        //h_ampl_rebin->Fit("gaus");
        //TFitResultPtr r = h_ampl_rebin->Fit("gaus");
        //r->Print();
        fout->Write();
        fout->Close();
        fclose(fp);
        return 0;
}
