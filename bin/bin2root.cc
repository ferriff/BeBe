#include "data_reader.h"

#include <iostream>

#include <getopt.h>

void help(const char * progr)
{
        std::cerr << "\n"
        "Usage: " << progr << " [options] [file(s)]\n"
        "\n"
        "where [options] can be:\n"
        "-h, --help                                print this help and exit\n"
        "-o <file>, --output=<file>                output root file, default: output.root\n"
        "-H <file>, --heat-file=<file>             trigger mode: input binary file, heat channel\n"
        "-L <file>, --light-file=<file>            trigger mode: input binary file, light channel\n"
        "-T <file>, --trigger-file=<file>          trigger mode: input text file, trigger information\n"
        //"-S <file>, --streamer-file=<file>         streamer mode: input binary file\n"
        "-s <nsamples>, --sample-number=<int>      number of samples for each pulse\n"
        "\n"
        "Converts binary files into pre-processed root files.\n"
        "The default is to read streamer-mode binary files [file(s)].\n"
        "If a trigger-mode binary has to be read, the three options -H, -L, -T have to be specified\n"
        "and possible arguments [file(s)] are ignored.\n"
        "\n";
        exit(0);
}


int main(int argc, char ** argv)
{
        static struct option long_options[] =
        {
                {"output"       , required_argument, 0, 'o'} , 
                {"heat-file"    , required_argument, 0, 'H'} , 
                {"light-file"   , required_argument, 0, 'L'} , 
                {"trigger-file" , required_argument, 0, 'T'} , 
                //{"streamer-file", required_argument, 0, 'S'} , 
                {"help"         , no_argument,       0, 'h'},
                {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;
        int c = 0, nsamples = 0;
        std::string ofile = "output.root", fheat, flight, ftrigger, fstreamer;
        int cnt = 0;
        while ((c = getopt_long(argc, argv, "o:H:L:T:S:s:h", long_options, &option_index)) != -1)
        {
                switch(c) {
                        case 'h':
                                help(argv[0]);
                                break;
                        case 'o':
                                ofile = optarg;
                                break;
                        case 'H':
                                fheat = optarg;
                                ++cnt;
                                break;
                        case 'L':
                                flight = optarg;
                                ++cnt;
                                break;
                        case 'T':
                                ftrigger = optarg;
                                ++cnt;
                                break;
                        case 's':
                                nsamples = atoi(optarg);
                                break;
                        default:
                                help(argv[0]);
                                break;
                }
        }
        bb::data_reader dr(ofile.c_str());
        if (nsamples) dr.set_n_samples(nsamples);
        if (cnt) {
                // trigger-mode binary file
                if (cnt != 3) {
                        std::cerr << "Error: Either all or none of -H, -T, -L must be specified.\n";
                        exit(-1);
                } else {
                        dr.read_trigger_mode_file(fheat.c_str(), flight.c_str(), ftrigger.c_str());
                }
        } else {
                // streamer-mode binary file
                for (int i = optind; i < argc; ++i) {
                        dr.read_streamer_mode_file(argv[i]);
                }
        }
        return 0;
}
