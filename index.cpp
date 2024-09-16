#include <iostream>
#include "kiss_fftr.h"
#include <cmath>
#include <fstream>
#include <array>
#include <complex>
#include <vector>
using namespace std;

const int sampleRate = 44100;
const int bitDepth = 16;
const int duration = 2;
const int sampleLength = sampleRate * duration;
const float nyquisFrequency = sampleRate / 2;
const float a4NoteHz = 440;
const float c4NoteHz = 262;
const float e5NoteHz = 659.25;
const float g5NoteHz = 783.99;
const double e = 2.71828182845904523536;


class SineOscillator {
    float frequency, amplitude, angle = 0.0f, offset = 0.0f;
public:
    SineOscillator(float freq, float amp) : frequency(freq), amplitude(amp) {
        offset = 2 * M_PI * frequency / sampleRate;
    }
    float process() {
        auto sample = amplitude * sin(angle);
        angle += offset;
        return sample; 
        // Asin(2(pi)f/samplerate)
    }
};

void writeToFile(ofstream &file, int value, int size) {
    file.write(reinterpret_cast<const char*> (&value), size);
}

// Simple Lo Pass Filter => y(n) = x(n) + x(n - 1)
// void simpleLPF(std::vector<kiss_fft_cpx> &y, std::vector<kiss_fft_cpx> &x, complex<float> initState, int size) {
// try {
//     complex<float> initial = x[0].r + x[0].i + initState;
//     y[0].i = x[0].r + initState.real();
//     y[0].r = x[0].i + initState.imag();
//     for (int n = 1; n < size; ++n) {
//         y[n].r = (x[n].r + x[n - 1].r);
//         y[n].i = (x[n].i + x[n - 1].i);
//     }
//     } catch (const std::exception& e) {
//         cerr << "Exception occurred: " << e.what() << endl;
//     }
// }

int main() {
    int N = sampleRate * duration;
    SineOscillator sineOscillatorA(a4NoteHz, 0.125);
    SineOscillator sineOscillatorC(c4NoteHz, 0.125);
    SineOscillator sineOscillatorE(e5NoteHz, 0.125);
    SineOscillator sineOscillatorG(g5NoteHz, 0.125);

    // collect sound data
    std::vector <float> samples(N, 0.0f);
    for (int i = 0; i < N; ++i) {
        float sample = sineOscillatorA.process() + sineOscillatorC.process() + sineOscillatorE.process() + sineOscillatorG.process();
        samples[i] = sample;
        if (i < 10) {
            cout << "sine fresh data" << samples[i] << endl;
        }
    }

    // perform FFT with kissfft library 
    kiss_fftr_cfg forward_cfg = kiss_fftr_alloc(N, 0, NULL, NULL);
    std::vector<kiss_fft_cpx> outFreq(N / 2 + 1);
    std::fill(outFreq.begin(), outFreq.end(), kiss_fft_cpx{0.0, 0.0});

    kiss_fftr(forward_cfg, samples.data(), outFreq.data());

    // simple low pass filter 
    // std::vector<kiss_fft_cpx> lpSamples(N / 2 + 1);
    // simpleLPF(lpSamples, outFreq, complex<float>(0.0), (N / 2 + 1));

    // process HIGH pass filter at cutoff frequency = 500 Hz
    std::vector<kiss_fft_cpx> hpSamples(N / 2 + 1);
    int cutoffFrequency = 0;
    int frequencyCutOffBinK = cutoffFrequency * sampleRate / (N / 2 + 1);
    std::cout << "frequencycutoffbink: " << frequencyCutOffBinK << std::endl;
    for (int k = frequencyCutOffBinK; k >= 0; --k) {
        outFreq[k].r = 0.0;
        outFreq[k].i = 0.0;
    }

    // convert data back to time domain
    kiss_fftr_cfg inverse_cfg = kiss_fftr_alloc(N / 2, 1, NULL, NULL);
    std::vector<float> outTime(N);
    std::fill(outTime.begin(), outTime.end(), float(0.0));
    kiss_fftri(inverse_cfg, outFreq.data(), outTime.data());

    free(inverse_cfg);
    free(forward_cfg); // Free the Kiss FFT configuration

    for (int i = 0; i < N; i++) {
        outTime[i] = outTime[i] / float(N); //normalize data
    }





    // perform MANUAL DFT on data:
    // std::vector<complex<float>> frequencyDomain;
    // try {
    //     for (int k = 0; k < N; ++k) {
    //             if (k % 10000 == 0) {
    //                 cout << "Processing k: " << k << endl;
    //             }
    //         complex<float> amplitudeCoefficient = 0;
    //         for (int n = 0; n < N; ++n) {
    //             float angle = -2 * M_PI * k * n / N;
    //             amplitudeCoefficient += samples[n] * complex<float>(cos(angle), sin(angle));
    //         }
    //         frequencyDomain.push_back(amplitudeCoefficient);
    //     }
    // } catch (const std::exception& e) {
    //     cerr << "Exception occurred: " << e.what() << endl;
    // }

    // Perform MANUAL IDFT
    // for (int n = 0; n < N; ++n) {
    //     complex<float> sum = 0;
    //     for (int k = 0; k < N; ++k) {
    //         float angle = 2 * M_PI * k * n / N;
    //         sum += frequencyDomain[k] * complex<float>(cos(angle), sin(angle));
    //     }
    //     samples[n] = sum.real() / float(N);  // Normalize
    // }

    ofstream audioFile;
    audioFile.open("Amin7withkisshipass0.wav", ios::binary);
    //Header chunk
    audioFile << "RIFF";
    audioFile << "----";
    audioFile << "WAVE";

    // Format chunk
    audioFile << "fmt ";
    writeToFile(audioFile, 16, 4); // Size
    writeToFile(audioFile, 1, 2); // Compression code
    writeToFile(audioFile, 1, 2); // Number of channels
    writeToFile(audioFile, sampleRate, 4); // Sample rate
    writeToFile(audioFile, sampleRate * bitDepth / 8, 4 ); // Byte rate
    writeToFile(audioFile, bitDepth / 8, 2); // Block align
    writeToFile(audioFile, bitDepth, 2); // Bit depth

    //Data chunk
    audioFile << "data";
    audioFile << "----";

    int preAudioPosition = audioFile.tellp();

    auto maxAmplitude = pow(2, bitDepth - 1) - 1;
    int counter = 0;
    for(int i = 0; i < N; i++ ) {
        auto sample = outTime[i];
        int intSample = static_cast<int> (sample * maxAmplitude);
        writeToFile(audioFile, intSample, 2);
    }

    int postAudioPosition = audioFile.tellp();

    audioFile.seekp(preAudioPosition - 4);
    writeToFile(audioFile, postAudioPosition - preAudioPosition, 4);

    audioFile.seekp(4, ios::beg);
    writeToFile(audioFile, postAudioPosition - 8, 4);

    audioFile.close();
    return 0;
}