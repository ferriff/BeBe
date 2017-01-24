#ifndef PULSE_H
#define PULSE_H
/*
 *        Class: Pulse
 *  Description: handling of sample and basic properties
 *       Author: Federico Ferri, CEA/Saclay
 */
#include "Types.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>

class Pulse
{
        public:
                Pulse(bbsize nsamples, const daqint * data = 0);

                ~Pulse();

                bbreal * data() { return _data; }
                //void print(FILE * fd);
                void printData(std::ostream & os);
                void inspect(std::ostream & os);

                void setData(const daqint * data);
                //void copyData(daqint * data) { }

                void preProcess(bbsize ped_samples);

                bbsize nSamples() { return _nsamples; };

                // pulse analysis
                bbreal average(bbsize start, bbsize size);
                bbreal rms(bbsize start, bbsize size);
                std::pair<bbreal, bbreal> maximum(bbsize start, bbsize size);
                std::pair<bbreal, bbreal> maximum_fitted(bbsize start, bbsize size);
                bbreal decay_time(bbsize imax, bbreal fraction);
                bbreal rise_time(bbsize imax, bbreal fraction);
                bbreal rise_time_interpolated(bbsize imax, bbreal fraction, bbreal amplitude, bbreal tmax);
                bbreal decay_time_interpolated(bbsize imax, bbreal fraction, bbreal amplitude, bbreal tmax);

        private:
                bbsize _nsamples;
                bbreal * _data;
                bbreal  * _fata;
};

#endif
