# BeBe: analysis framework for double beta decay bolometers

### Installation instructions

#### Requirements:
   * `ROOT` package installed

#### Steps:
```
git clone git@github.com:ferriff/BeBe.git # or, if no github account `git clone https://github.com/ferriff/BeBe.git'
cd BeBe
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
     ```
size_t fignal_rise_time(const float data[], size_t imax, float fraction)
{
        float m = data[imax];
        for (size_t i = imax; i >= 0; --i) {
                if (data[i] / m < fraction) return imax - i;
        }
        return 0;
}
     ```
