#include "pulse.h"

#include <cmath>
#include <fstream>
#include <iostream>

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
        bb::daqint_t * data;
        Double_t ts;
} Event;

Event  _e;

void set_branches(TTree * t)
{
        t->SetBranchAddress("detid", &_e.detid);
        // handle variable size array
        // allocate maximum dimension (cf. TTreePlayer::MakeClass code)
        TLeaf * l = t->GetLeaf("nsamples", "nsamples");
        _e.data = (bb::daqint_t *)calloc(l->GetMaximum(), sizeof(bb::daqint_t));
        t->SetBranchAddress("raw_pulse", _e.data);
        t->SetBranchAddress("nsamples", &_e.nsamples);
        t->SetBranchAddress("detid", &_e.detid);
        t->SetBranchAddress("time", &_e.ts);
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
        TH1F     * h_ampl_res               = new TH1F("h_ampl_res", "h_ampl_res", 1000, -10., 10.);
        TH1F     * h_ped                    = new TH1F("h_ped", "h_ped", 2<<13, 0., 2<<17);
        TH1F     * h_ped_raw                = new TH1F("h_ped_raw", "h_ped_raw", 1000, -1e5, 1e5);
        TH1F     * h_ped_rms                = new TH1F("h_ped_rms", "h_ped_rms", 100, 0., 100.);
        TH1F     * h_dt                     = new TH1F("h_dt", "h_dt", 5000, 0., 5000000.);
        TH1F     * h_time_s1                = new TH1F("h_time_s1", "h_time_s1", 20000, 0., 50.);
        TGraph   * g_decay_vs_rise          = new TGraph();
        g_decay_vs_rise->SetNameTitle("g_decay_vs_rise", "g_decay_vs_rise");
        gDirectory->Add(g_decay_vs_rise);
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
        TGraph * g_rise_vs_ampl        = new TGraph();
        g_rise_vs_ampl->SetNameTitle("g_rise_vs_ampl", "g_rise_vs_ampl");
        gDirectory->Add(g_rise_vs_ampl);
        TGraph * g_decay_vs_ampl        = new TGraph();
        g_decay_vs_ampl->SetNameTitle("g_decay_vs_ampl", "g_decay_vs_ampl");
        gDirectory->Add(g_decay_vs_ampl);
        TFile * fin = TFile::Open("lsm_signals.root");
        TTree * t = (TTree*)fin->Get("ntu");
        fout->cd();
        set_branches(t);
        Long64_t nentries = t->GetEntries();
        Long64_t gcnt = 0, gcnt2 = 0, gcnta = 0, gcnt_rd = 0;
        UInt_t otmaxdaq = 1;
        FILE * fp = fopen("pulses.dat", "w");
        size_t ipulse = 0, totpulse = 1;
        std::ofstream ofs, pfs;
        ofs.open ("pappa.dat", std::ofstream::out);
        pfs.open ("pippa.dat", std::ofstream::out);
        for (Long64_t ien = 0; ien < nentries; ++ien) {
                t->GetEntry(ien);
                // select only one channel
                if (_e.detid != 5) continue;
                bb::pulse p(_e.nsamples);
                p.set_data(_e.data);
                /*
                h_ampl_raw->Fill(_e.ampl);
                */
                // signal analysis
                h_ped_raw->Fill(p.average(0, 50));
                p.pre_process(100);
                //p.print_data(std::cerr);
                float ped = p.average(0, 50);
                float ped_rms = p.rms(0, 50);
                h_ped->Fill(ped);
                h_ped_rms->Fill(ped_rms);
                // detailed check of pulses if conditions applies
                /*
                if (ped_rms > 10) {
                        p.inspect(ofs);
                        ofs << "# pulse number:" << ipulse << "\n";
                        ofs << "\n\n";
                        ofs.flush();
                        // pause until Enter key pressed
                        std::cerr << "pulse " << ipulse << " of " << totpulse << " dumped, press [Enter] to continue...\n";
                        //getchar(); // uncomment if you want to pause
                        ++ipulse;
                        continue; // do not perform the rest of the analysis
                }
                */
                auto res = p.maximum(100, p.n_samples());
                size_t iM = res.first;
                float M = res.second;
                res = p.maximum_fitted(100, p.n_samples());
                float fiM = res.first;
                float fM = res.second;
                h_ampl->Fill(fM);
                h_ampl_res->Fill((fM - M) / fM);
                //std::cerr << "--> " << iM << " " << M << " " << fiM << " " << fM << "\n";
                float ft_daq = _e.ts + fiM; // in seconds
                float trise = p.rise_time_interpolated(iM, 0.05, fM, fiM) - p.rise_time_interpolated(iM, 0.95, fM, fiM);
                float tdecay = p.decay_time_interpolated(iM, 0.20, fM, fiM) - p.decay_time_interpolated(iM, 0.95, fM, fiM);
                if(!(std::isnan(trise) || std::isnan(tdecay)) && fM > 0) {
                        //fprintf(stderr, "times %f %f %f\n", fM, trise, tdecay);
                        if (trise > 30 && trise < 31 && fM > 1e4 && fM < 1.1e4) {
                                p.inspect(ofs);
                                ofs << "# pulse number:" << ipulse << "\n";
                                ofs << "\n\n";
                                ofs.flush();
                        }
                        if (trise < 33 && trise > 32 && fM > 1e4 && fM < 1.1e4) {
                                p.inspect(pfs);
                                pfs << "# pulse number:" << ipulse << "\n";
                                pfs << "\n\n";
                                pfs.flush();
                        }
                        g_rise_vs_ampl->SetPoint(gcnt_rd, fM, trise);
                        g_decay_vs_ampl->SetPoint(gcnt_rd, fM, tdecay);
                        g_decay_vs_rise->SetPoint(gcnt_rd++, trise, tdecay);
                }
                /*
                fprintf(stdout, "%f %f %lu %u %u %lu %lu %f %f %f %f %f %f\n",
                        _e.ampl, M, iM, _e.t_max, _e.t_max_daq, 
                        signal_rise_time(fata, iM, 0.05), signal_decay_time(fata, iM, 0.20),
                        trise,
                        tdecay, 
                        fM, fiM, ped, ped_rms);
                */
                float adc2keV = 4783. / 101279.;
                h_ampl_keV->Fill(fM * adc2keV);
                h_ampl_drift_corr->Fill(fM - 1.788 * ft_daq / 3600.);
                //if (otmaxdaq && _e.ampl > 101100 && _e.ampl < 101500) {
                if (otmaxdaq && fM > 46400 && fM < 46500) {
                        h_dt->Fill(_e.ts - otmaxdaq);
                        h_ampl_heater->Fill(fM);
                        h_ampl_heater_drift_corr->Fill(fM - 1.788 * ft_daq / 3600.);
                        otmaxdaq = _e.ts;
                }
                if (fM > 46400 && fM < 46500) {
                        g_ampl_vs_t->SetPoint(gcnt++, ft_daq / 3600., fM);
                } else if (fM > 100.8e+3 && fM < 101.6e+3) {
                        g_ampl2_vs_t->SetPoint(gcnt2++, ft_daq / 3600., fM);
                }
                g_ampl_all_vs_t->SetPoint(gcnta++, ft_daq / 3600., fM);
                if (fM > 1e5) h_time_s1->Fill(trise - 0.01 * fM / 1000.);
                //else if (_e.ampl > 101.2e+3 && _e.ampl < 101.4e+3) g_ampl2_vs_t->SetPoint(gcnt2++, ft_daq / 2000. / 3600., fM);
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
