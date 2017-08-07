#ifndef EVENT_H
#define EVENT_H
/*
 *        Class: event
 *  Description: single detector event description
 *       Author: Federico Ferri, CEA/Saclay
 */

#include "types.h"

#include "TFile.h"

namespace bb {

        class event {

                friend class tree_reader;

                public:

                        ~event();

                        daqint_t * pulse() const { return _data; }
                        uint_t     nsamples() const { return _nsamples; }
                        int32_t    detid() const { return _detid; }
                        Double_t   time()  const { return _time; }
                        //Long64_t * event_ids() const { Long64_t * dest = (Long64_t*)malloc(sizeof(_event_ids)); return (Long64_t*)memcpy(dest, _event_ids, sizeof(_event_ids) / sizeof(Long64_t)); }
                        Long64_t * event_ids() const { return _event_ids; }

                private:
                        void init_event_ids(size_t size);
                        void init_data(size_t size);
                        void set_data_size(size_t size);

                        int32_t _detid;
                        unsigned short _nsamples;
                        daqint_t * _data = 0;
                        Double_t _time;
                        Double_t _freq; 
                        Long64_t * _event_ids = 0;
        };
}

#endif
