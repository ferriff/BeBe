#include "tree_reader.h"
#include "types.h"

void print_pulse(bb::event & e)
{
        bb::daqint_t * pulse = e.pulse();
        bb::uint_t ns = e.nsamples();
        for (unsigned int i = 0; i < ns; ++i) {
                fprintf(stderr, "%u  %d  %d\n", i, pulse[i], e.detid());
        }
}

int main()
{
        //bb::tree_reader tr("pappa_derivative.root");
        bb::tree_reader tr("test2.root");
        tr.set_detector("LMO1t_chaleur");
        auto curr_det = tr.detector();
        fprintf(stderr, "_ct = %lu -> %s\n", tr.detector(), tr.detector_name().c_str());
        //tr.set_detector("CWO_chaleur");
        size_t cnt = 0;
        Long64_t * ids = (Long64_t*)malloc(sizeof(Long64_t) * tr.n_detectors());
        while (tr.next_event()) {
                auto t = tr.e().time();
                // copy to avoid that changing the detector points to a different set of ids
                memcpy(ids, tr.e().event_ids(), sizeof(Long64_t) * tr.n_detectors());
                fprintf(stderr, "--> %f   ", t);
                for (size_t i = 0; i < tr.n_detectors(); ++i) {
                        fprintf(stderr, " %2lu:%lld", i, ids[i]);
                }
                fprintf(stderr, "\n");
                //print_pulse(tr.e());
                for (size_t i = 0; i < tr.n_detectors(); ++i) {
                        tr.set_detector(i);
                        auto ret_prev = tr.read_event(ids[i]);
                        auto t_prev = tr.e().time();
                        if (ret_prev == 0) t_prev = -1;
                        auto ret_next = tr.read_event(ids[i] + 1);
                        auto t_next = tr.e().time();
                        if (ret_next == 0) t_next = -1;
                        fprintf(stderr, ":=: %2lu %15s  %+9.4f %8lld    %+9.4f    %2d %2d\n", i, tr.detector_name().c_str(), t_prev, ids[i], t_next, ret_prev, ret_next);
                }
                tr.set_detector(curr_det);
                ++cnt;
        }
        free(ids);
        fprintf(stderr, "end> %f\n", tr.e().time());
        return 0;
}
