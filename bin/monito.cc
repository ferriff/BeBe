#include "tree_reader.h"

#include "pulse.h"
#include "types.h"

#include "TGraph.h"
#include "TH1F.h"
#include "TProfile.h"

#include "HistoManager.h"

#define TIME_WINDOW 3000 // FIXME

HistoManager hm_;

void print_pulse(bb::event & e)
{
        bb::daqint_t * pulse = e.pulse();
        bb::uint_t ns = e.nsamples();
        for (unsigned int i = 0; i < ns; ++i) {
                fprintf(stderr, "%u  %d  %d\n", i, pulse[i], e.detid());
        }
}


void process_detector(bb::tree_reader & tr, const char * name)
{
        tr.set_detector(name);
        auto curr_det = tr.detector();
        fprintf(stderr, "analyzing detector %s %lu\n", tr.detector_name().c_str(), tr.detector());
        size_t cnt = 0;
        Long64_t * ids = (Long64_t*)malloc(sizeof(Long64_t) * tr.n_detectors());
        while (tr.next_event()) {

                auto & e = tr.e();
                bb::pulse p(e.nsamples(), e.pulse());

                float ped_avg = p.average(0, 200);
                float ped_rms = p.rms(0, 200);
                static std::string str_ped = std::string("pedestal_") + name;
                hm_.h<TH1F>("amplitude", str_ped.c_str())->Fill(ped_avg);
                hm_.h<TH1F>("pedestalrms", name)->Fill(ped_rms);
                // subtract the pedestal
                // (no additional filters applied for the time being)
                p.pre_process(200);

                // look for the maximum within a window of 512 samples around
                // the triggering one (at nsamples / 3)
                auto res = p.maximum(e.nsamples() / 3 - 256, e.nsamples() / 3 + 256);
                //float max_idx  = res.first;
                float max_ampl = res.second;
                hm_.h<TH1F>("amplitude", name)->Fill(max_ampl);

                // find the minimum to reject pulses with spikes
                res = p.minimum(0, TIME_WINDOW);
                float min_ampl = res.second;
                static std::string str_min = std::string("minimum_") + name;
                hm_.h<TH1F>("amplitude", str_min.c_str())->Fill(min_ampl);

                // average pulse
                bb::real_t * d = p.data();
                for (size_t i = 0; i < e.nsamples(); ++i) {
                        hm_.h<TProfile>("pulse", name)->Fill(i, d[i]);
                }

                // fitted amplitude
                res = p.maximum_fitted(100, p.n_samples());
                float max_idx_fit  = res.first;
                float max_ampl_fit = res.second;
                static std::string str_fit = std::string("fitted_") + name;
                hm_.h<TH1F>("amplitude", str_fit.c_str())->Fill(max_ampl_fit);
                static std::string str_time = std::string("maxtime_") + name;
                hm_.h<TH1F>("hpulse", str_time.c_str())->Fill((float)max_idx_fit);

                tr.set_detector(curr_det);
                ++cnt;
        }
        free(ids);
}


int main()
{
        // signal and noise average pulses
        hm_.addTemplate<TProfile>("pulse", new TProfile("p_pulse","p_pulse", TIME_WINDOW, 0., TIME_WINDOW));
        hm_.addTemplate<TH1F>("hpulse", new TH1F("h_pulse","h_pulse", TIME_WINDOW, 0., TIME_WINDOW));
        // amplitudes
        hm_.addTemplate<TH1F>("amplitude", new TH1F("h_amplitude", "h_amplitude", 40000, -20000, 20000));
        hm_.addTemplate<TH1F>("resolution", new TH1F("h_resolution", "h_resolution", 1000, -10., 10.));
        hm_.addTemplate<TH1F>("pedestalrms", new TH1F("h_pedestal_rms", "h_pedestal_rms", 100, 0., 100.));
        // for holding graphs
        auto g = new TGraph();
        g->SetNameTitle("g_graph", "g_graph");
        hm_.addTemplate<TGraph>("graph", std::move(g));
        hm_.print();

        bb::tree_reader tr("test2.root");
        for (auto & name : tr.detector_names()) {
                process_detector(tr, name.c_str());
        }
        hm_.save("histos.root");

        return 0;
}
