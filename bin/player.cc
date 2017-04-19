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
        bb::tree_reader tr("pappa_derivative.root");
        tr.set_detector(6);
        size_t cnt = 0;
        while (tr.next_event()) {
                fprintf(stderr, "--> %f\n", tr.e().time());
                print_pulse(tr.e());
                for (int i = 0; i < 12; ++i) {
                }
                break;
        }
        return 0;
}
