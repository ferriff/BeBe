#include "event.h"

bb::event::~event()
{
        if (_data) free(_data);
        if (_event_ids) free(_event_ids);
}


void bb::event::init_event_ids(size_t size)
{
        if (_event_ids) free(_event_ids);
        _event_ids = (Long64_t *)calloc(size, sizeof(Long64_t));
}


void bb::event::init_data(size_t size)
{
        if (_data) free(_data);
        set_data_size(size);
}


void bb::event::set_data_size(size_t size)
{
        if (sizeof(_data) / sizeof(daqint_t) >= size) return;
        _data = (daqint_t *)calloc(size, sizeof(daqint_t));
}
