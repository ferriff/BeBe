#include "Pulse.h"

#include <iostream>
#include <fstream>

#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLeaf.h"
#include "TProfile.h"
#include "TTree.h"

#define TIME_WINDOW 2048

typedef struct {
        int32_t detid;
        Int_t  nsamples;
        daqint * data;
} Event;

Event  _e;
//size_t _ipulse;

void set_branches(TTree * t)
{
        t->SetBranchAddress("detid", &_e.detid);
        // handle variable size array
        // allocate maximum dimension (cf. TTreePlayer::MakeClass code)
        TLeaf * l = t->GetLeaf("nsamples", "nsamples");
        _e.data = (daqint *)calloc(l->GetMaximum(), sizeof(daqint));
        t->SetBranchAddress("raw_pulse", _e.data);
        /////t->SetBranchAddress("amplitude", &_e.ampl);
        /////t->SetBranchAddress("baseline", &_e.bl);
        /////t->SetBranchAddress("t_max", &_e.t_max);
        /////t->SetBranchAddress("t_max_daq", &_e.t_max_daq);
        /////t->SetBranchAddress("t_start_daq", &_e.t_start_daq);
        /////t->SetBranchAddress("pu", &_e.pu);
        t->SetBranchAddress("nsamples", &_e.nsamples);
        //return br;
}


int main()
{
        TFile * fout = TFile::Open("histos.root", "recreate");
        TH1F     * h_ampl_raw               = new TH1F("h_ampl_raw", "h_ampl_raw", 2<<17, 0., 2<<17);
        TH1F     * h_ampl                   = (TH1F*)h_ampl_raw->Clone("h_ampl");
        TH1F     * h_ampl_drift_corr        = (TH1F*)h_ampl_raw->Clone("h_ampl_drift_corr");
        TH1F     * h_ampl_heater            = (TH1F*)h_ampl_raw->Clone("h_ampl_heater");
        TH1F     * h_ampl_heater_drift_corr = (TH1F*)h_ampl_raw->Clone("h_ampl_heater_drift_corr");
        TH1F     * h_ampl_keV               = new TH1F("h_ampl_keV", "h_ampl_keV", 4000, 0., 10000.);
        TH1F     * h_ped                    = new TH1F("h_ped", "h_ped", 2<<13, 0., 2<<17);
        TH1F     * h_ped_raw                = new TH1F("h_ped_raw", "h_ped_raw", 1000, -1e5, 1e5);
        TH1F     * h_ped_rms                = new TH1F("h_ped_rms", "h_ped_rms", 100, 0., 100.);
        TH1F     * h_dt                     = new TH1F("h_dt", "h_dt", 5000, 0., 5000000.);
        TH1F     * h_time_s1                = new TH1F("h_time_s1", "h_time_s1", 20000, 0., 50.);
        TProfile * p_average_signal         = new TProfile("p_average_signal", "p_average_signal", 5. * TIME_WINDOW, 0., TIME_WINDOW);
        TGraph * g_baseline_vs_t            = new TGraph();
        g_baseline_vs_t->SetNameTitle("g_baseline_vs_t", "g_baseline_vs_t");
        gDirectory->Add(g_baseline_vs_t);
        TGraph * g_ampl_vs_t            = new TGraph();
        g_ampl_vs_t->SetNameTitle("g_ampl_vs_t", "g_ampl_vs_t");
        gDirectory->Add(g_ampl_vs_t);
        TGraph * g_ampl2_vs_t = new TGraph();
        g_ampl2_vs_t->SetNameTitle("g_ampl2_vs_t", "g_ampl2_vs_t");
        gDirectory->Add(g_ampl2_vs_t);
        TGraph * g_ampl_all_vs_t        = new TGraph();
        g_ampl_all_vs_t->SetNameTitle("g_ampl_all_vs_t", "g_ampl_all_vs_t");
        gDirectory->Add(g_ampl_all_vs_t);
        TFile * fin = TFile::Open("lsm_signals.root");
        TTree * t = (TTree*)fin->Get("ntu");
        fout->cd();
        set_branches(t);
        Long64_t nentries = t->GetEntries();
        Long64_t gcnt = 0, gcnt2 = 0, gcnta = 0;
        UInt_t otmaxdaq = 1;
        float fata[TIME_WINDOW]; // FIXME
        FILE * fp = fopen("pulses.dat", "w");
        size_t ipulse = 1, totpulse = 1;
        std::ofstream ofs;
        ofs.open ("pappa.dat", std::ofstream::out);
        for (Long64_t ien = 0; ien < nentries; ++ien) {
                t->GetEntry(ien);
                // select only one channel
                if (_e.detid != 5) continue;
                Pulse p(_e.nsamples);
                p.setData(_e.data);
                /*
                h_ampl_raw->Fill(_e.ampl);
                */
                // signal analysis
                h_ped_raw->Fill(p.average(0, 50));
                p.preProcess(100);
                //p.printData(std::cerr);
                float ped = p.average(0, 50);
                float ped_rms = p.rms(0, 50);
                h_ped->Fill(ped);
                h_ped_rms->Fill(ped_rms);
                // detailed check of pulses if conditions applies
                if (ped_rms > 10) {
                        p.inspect(ofs);
                        ofs << "# pulse number:" << ipulse << "\n";
                        ofs << "\n\n";
                        ofs.flush();
                        // pause until Enter key pressed
                        std::cerr << "pulse " << ipulse << " of " << totpulse << " dumped, press [Enter] to continue...\n";
                        //getchar();
                        ++ipulse;
                }
                auto res = p.maximum(100, p.nSamples());
                size_t iM = res.first;
                float M = res.second;
                res = p.maximum_fitted(100, p.nSamples());
                float fiM = res.first;
                float fM = res.second;
                h_ampl->Fill(fM);
                /*
                float ft_daq = _e.t_max_daq - (_e.t_max - iM);
                float trise = signal_rise_time_interpolated(fata, iM, 0.05, fM, fiM) - signal_rise_time_interpolated(fata, iM, 0.95, fM, fiM);
                float tdecay = signal_decay_time_interpolated(fata, iM, 0.20, fM, fiM) - signal_decay_time_interpolated(fata, iM, 0.95, fM, fiM);
                fprintf(stdout, "%f %f %lu %u %u %lu %lu %f %f %f %f %f %f\n",
                        _e.ampl, M, iM, _e.t_max, _e.t_max_daq, 
                        signal_rise_time(fata, iM, 0.05), signal_decay_time(fata, iM, 0.20),
                        trise,
                        tdecay, 
                        fM, fiM, ped, ped_rms);
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
                */
                ++totpulse;
        }
        ofs.close();
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
