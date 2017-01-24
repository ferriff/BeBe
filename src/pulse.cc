#include "Pulse.h"

#include "TF1.h"
#include "TGraph.h"

#include <algorithm>
#include <limits>

Pulse::Pulse(bbsize nsamples, const daqint * data) :
        _nsamples(nsamples)
{
        _data = (bbreal *)calloc(_nsamples, sizeof(bbreal));
        if (data) setData(data);
}

Pulse::~Pulse()
{
        //if (_data) free(_data);
}

void Pulse::setData(const daqint * data)
{
        for (bbsize i = 0; i < _nsamples; ++i) {
                _data[i] = (bbreal)data[i];
        }
}


//void Pulse::print(FILE * fd)
//{
//        for (bbsize i = 0; i < _nsamples; ++i) {
//                fprintf(fd, "%u %d\n", i, _data[i]); // FIXME
//        }
//}


void Pulse::printData(std::ostream & os)
{
        for (bbsize i = 0; i < _nsamples; ++i) {
                os << i << " " << _data[i] << "\n";
        }
}


void Pulse::inspect(std::ostream & os)
{
        printData(os);
        os << "#        average: " << average(0, _nsamples) << "\n";
        os << "#  average[0:50]: " << average(0, std::min((bbsize)50, _nsamples)) << "\n";
        os << "# average[0:100]: " << average(0, std::min((bbsize)100, _nsamples)) << "\n";
        os << "#            rms: " << rms(0, _nsamples) << "\n";
        os << "#      rms[0:50]: " << rms(0, std::min((bbsize)50, _nsamples)) << "\n";
        auto p = maximum(0, _nsamples);
        os << "#            max: " << p.second << " at " << p.first << "\n";
        p = maximum_fitted(0, _nsamples);
        os << "#     fitted max: " << p.second << " at " << p.first << "\n";
}


void Pulse::preProcess(bbsize ped_samples)
{
        bbreal ped = average(0, ped_samples);
        for (bbsize i = 0; i < _nsamples; ++i) _data[i] = _data[i] - ped;
}


bbreal Pulse::average(bbsize start, bbsize size)
{
        bbreal m = 0;
        for (bbsize i = start; i < size; ++i) {
                m += _data[i];
        }
        return (bbreal)m / (bbreal)size;
}


bbreal Pulse::rms(bbsize start, bbsize size)
{
        bbreal m = 0, mm = 0;
        for (bbsize i = start; i < size; ++i) {
                m  += _data[i];
                mm += _data[i] * _data[i];
        }
        m /= (bbreal)(size - start);
        return sqrt(mm / (bbreal)(size - start) - m * m);
}


std::pair<bbreal, bbreal> Pulse::maximum(bbsize start, bbsize size)
{
        bbreal m = -std::numeric_limits<bbreal>::max();
        //std::cerr << "--> m = " << m << "\n";
        //fprintf(stderr, "--> m = %lf\n", m);
        bbreal im = 0;
        for (bbsize i = start; i < size; ++i) {
                //fprintf(stderr, "--> %d %f %f\n", i, m, _data[i]);
                if (_data[i] > m) {
                        m = _data[i];
                        im = i;
                }
                //getchar();
        }
        return std::make_pair(im, m);
}


std::pair<bbreal, bbreal> Pulse::maximum_fitted(bbsize start, bbsize size)
{
        bbreal m = -std::numeric_limits<bbreal>::max();
        bbsize im = 0;
        for (bbsize i = start; i < size; ++i) {
                if (_data[i] > m) {
                        m = _data[i];
                        im = i;
                }
        }
        //char tmp[16];
        //sprintf(tmp, "pulse_%06lu", _ipulse);
        TGraph * g = new TGraph();
        //g->SetNameTitle(tmp, tmp);
        int npoints = 10;
        for (int i = 0; i < npoints; ++i) g->SetPoint(i, (bbreal)im - npoints / 2 + i, (bbreal)_data[im - npoints / 2 + i]);
        g->Fit("pol2", "Q");
        TF1 * f = g->GetFunction("pol2");
        bbreal p0 = f->GetParameter(0);
        bbreal p1 = f->GetParameter(1);
        bbreal p2 = f->GetParameter(2);
        delete g;
        // return tmax, max
        return std::make_pair(-0.5 * p1 / p2, p0 - 0.25 * p1 * p1 / p2);
}


bbreal Pulse::decay_time(bbsize imax, bbreal fraction)
{
        bbreal m = _data[imax];
        for (bbsize i = imax; i < _nsamples; ++i) {
                if (_data[i] / m < fraction) return i - imax;
        }
        return 0;
}


bbreal Pulse::rise_time(bbsize imax, bbreal fraction)
{
        bbreal m = _data[imax];
        for (bbsize i = imax; i >= 0; --i) {
                if (_data[i] / m < fraction) return imax - i;
        }
        return 0;
}


bbreal Pulse::rise_time_interpolated(bbsize imax, bbreal fraction, bbreal amplitude, bbreal tmax)
{
        bbreal m = _data[imax];
        if (amplitude) m = amplitude;
        for (bbsize i = imax; i >= 0; --i) {
                if (_data[i] / m < fraction) {
                        //return imax - i + (m * fraction - _data[i]) / (_data[i + 1] - _data[i]) * 1.;
                        return tmax - (i + 1 + (_data[i + 1] - m * fraction) / (_data[i + 1] - _data[i]) * 1.);
                }
        }
        return -1;
}


bbreal Pulse::decay_time_interpolated(bbsize imax, bbreal fraction, bbreal amplitude, bbreal tmax)
{
        bbreal m = _data[imax];
        if (amplitude) m = amplitude;
        for (bbsize i = imax; i < _nsamples; ++i) {
                if (_data[i] / m < fraction) {
                        return i - (m * fraction - _data[i]) / (_data[i - 1] - _data[i]) - tmax;
                }
        }
        return -1;
}
