#include "DataReader.h"

int main(int argc, char ** argv)
{
        DataReader dr("pappa.root");
        for (int i = 1; i < argc; ++i) {
                dr.readStreamerModeFile(argv[i]);
        }
        return 0;
}
