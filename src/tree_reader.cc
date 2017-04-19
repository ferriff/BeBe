#include "tree_reader.h"

bb::tree_reader::tree_reader(const char * filename)
{
        _f = TFile::Open(filename);
        if (!_f || _f->IsZombie()) {
                fprintf(stderr, "[bb::tree_reader::tree_reader] problems in opening input file `%s', aborting.\n", filename);
                exit(-1);
        }
        TTree * t_md = (TTree*)_f->Get("MetaData/meta");
        t_md->SetBranchAddress("tree_order", &_t_names);
        t_md->GetEntry(0);
        // assign tree in the same order as they were read from the raw data
        // so that event_ids are set up correctly
        for (auto & n : *_t_names) {
                _t.push_back((TTree *)_f->Get((std::string("Data/") + n).c_str()));
                _t_entry.push_back(0);
        }
        _evt.init_event_ids(_t.size());
        init_trees();
}


void bb::tree_reader::init_trees()
{
        for (auto & t : _t) {
                init_branches(t);
        }
}


void bb::tree_reader::init_branches(TTree * t)
{
        t->SetBranchAddress("nsamples", &_evt._nsamples);
        t->GetEntry(0);
        _evt.init_data(_evt._nsamples * 5);
        t->SetBranchAddress("raw_pulse", _evt._data);
        t->SetBranchAddress("detid", &_evt._detid);
        t->SetBranchAddress("time", &_evt._time);
        t->SetBranchAddress("freq", &_evt._freq);
        char tmp[32];
        sprintf(tmp, "event_ids[%lu]", _t.size());
        t->SetBranchAddress(tmp, _evt._event_ids);
}


bool bb::tree_reader::next_event()
{
        fprintf(stderr, "--> %d %u\n", _ct, _t_entry[_ct]);
        return read_event(_t_entry[_ct]++);
}


bool bb::tree_reader::prev_event()
{
        return read_event(_t_entry[_ct]--);
}


bool bb::tree_reader::read_event(ULong64_t entry)
{
        fprintf(stderr, "--> %p %u\n", _t[_ct], entry);
        return _t[_ct]->GetEntry(entry);
}


void bb::tree_reader::set_detector(const char * dname)
{
        auto it = std::find(_t_names->begin(), _t_names->end(), dname);
        if (it != _t_names->end()) _ct = it - _t_names->begin();
        else {
                fprintf(stderr, "[bb::tree_reader::set_detector] Warning: detector `%s' not found, reset to detector 0.\n", dname);
                _ct = 0;
        }
}


void bb::tree_reader::set_detector(size_t idet)
{
        _ct = idet;
}
