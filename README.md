# bebe

Analysis framework for double beta decay bolometers

### Installation instructions

##### Requirements:
   * `CMake` package (>= 3.4.3), see https://cmake.org/
   * `ROOT` package, see http://root.cern.ch

##### Steps:

```
git clone git@github.com:ferriff/bebe.git
# or, if no github account `git clone https://github.com/ferriff/bebe.git'
cd bebe
cmake .
make
```

### Project structure

   * `include/`: header files for libraries
   * `src/`: class implementation files for libraries
   * `bin/`: analysis codes

Whenever a new file is added, remember to update `src/CMakeLists.txt`
and/or `bin/CMakeLists.txt` accordingly. Follow the example of what
already there.

### Minimal coding rules:
   * `.h` suffix for headers, `.cc` for implementation and analysis code
   * 8 space indentation (no TABs)
   * example snapshot
```C++
size_t fignal_rise_time(const float data[], size_t imax, float fraction)
{
        float m = data[imax];
        for (size_t i = imax; i >= 0; --i) {
                if (data[i] / m < fraction) return imax - i;
        }
        return 0;
}
```
