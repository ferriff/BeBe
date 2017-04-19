#ifndef TREE_READER_H
#define TREE_READER_H
/*
 *        Class: tree_reader
 *  Description: handling of tree files with multiple detectors
 *       Author: Federico Ferri, CEA/Saclay
 */
#include "event.h"

#include <vector>

#include "TFile.h"
#include "TTree.h"

namespace bb {

        class tree_reader
        {
                public:
                        tree_reader(const char * filename);

                        ~tree_reader() {};

                        event & e() { return _evt; }

                        void set_detector(const char * dname);
                        void set_detector(size_t idet);

                        bool next_event();
                        bool prev_event();
                        bool read_event(ULong64_t entry);

                private:
                        void init_branches(TTree * t);
                        void init_trees();

                        bb::event _evt;
                        TFile * _f;
                        std::vector<TTree *> _t;
                        std::vector<ULong64_t> _t_entry;
                        std::vector<std::string> * _t_names = 0;
                        int _ct = 0; // current tree index
        };
}

#endif
