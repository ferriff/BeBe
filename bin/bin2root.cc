#include "data_reader.h"

int main(int argc, char ** argv)
{
        bb::data_reader dr("pappa.root");
        for (int i = 1; i < argc; ++i) {
                dr.read_streamer_mode_file(argv[i]);
        }
        return 0;
}
